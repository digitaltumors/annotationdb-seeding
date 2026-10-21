import csv
import json
import time

import requests


ANNOTATIONDB_URL = "http://127.0.0.1:8000/compound/many"
BATCH_SIZE = 200
NUMBER_OF_BATCHES = 16


def extract_unique_cids(atc_data):
    cids = set()

    for annotation in atc_data.get("Annotations", []):
        linked_records = annotation.get("LinkedRecords", {})

        for cid in linked_records.get("CID", []):
            cids.add(int(cid))

    return sorted(cids)


def build_cid_to_atc(atc_data):
    """Build a mapping of CID -> semicolon-separated ATC codes."""
    cid_to_codes = {}

    for annotation in atc_data.get("Annotations", []):
        atc_code = annotation.get("SourceID", "")
        if not atc_code:
            continue

        linked_records = annotation.get("LinkedRecords", {})

        for cid in linked_records.get("CID", []):
            cid = int(cid)
            if cid not in cid_to_codes:
                cid_to_codes[cid] = []
            if atc_code not in cid_to_codes[cid]:
                cid_to_codes[cid].append(atc_code)

    # Join multiple codes with semicolons
    return {
        cid: ";".join(sorted(codes))
        for cid, codes in cid_to_codes.items()
    }


def extract_returned_cids(response_json):
    if isinstance(response_json, list):
        compounds = response_json
    elif isinstance(response_json, dict):
        compounds = (
            response_json.get("compounds")
            or response_json.get("results")
            or response_json.get("data")
            or []
        )
    else:
        compounds = []

    returned_cids = set()

    for compound in compounds:
        if not isinstance(compound, dict):
            continue

        cid = compound.get("cid")

        if cid is not None:
            returned_cids.add(int(cid))

    return returned_cids


def request_batch(cid_batch, maximum_attempts=5):
    params = []

    for cid in cid_batch:
        params.append(("compound", str(cid)))

    params.extend([
        ("format", "json"),
        ("bioassay", "true"),
        ("mechanism", "true"),
        ("toxicity", "true"),
        ("golden_bioassay", "true"),
    ])

    for attempt in range(1, maximum_attempts + 1):
        try:
            response = requests.get(
                ANNOTATIONDB_URL,
                params=params,
                timeout=180,
            )

            response.raise_for_status()
            return response

        except requests.RequestException as error:
            print(
                f"Attempt {attempt}/{maximum_attempts} failed: {error}"
            )

            if attempt == maximum_attempts:
                raise

            time.sleep(attempt * 3)


def main():
    with open("atc_data.json", "r", encoding="utf-8") as file:
        atc_data = json.load(file)

    all_cids = extract_unique_cids(atc_data)
    cid_to_atc = build_cid_to_atc(atc_data)

    requested_cids = all_cids[
        :BATCH_SIZE * NUMBER_OF_BATCHES
    ]

    batch_summaries = []
    all_returned_cids = set()

    print(f"Unique CIDs in input: {len(all_cids)}")
    print(f"CIDs being requested: {len(requested_cids)}")

    total_batch = []
    for batch_index in range(NUMBER_OF_BATCHES):
        batch_number = batch_index + 1
        start = batch_index * BATCH_SIZE
        end = start + BATCH_SIZE

        cid_batch = requested_cids[start:end]

        response = request_batch(cid_batch)
        response_json = response.json()

        returned_cids = extract_returned_cids(response_json)
        returned_requested_cids = set(cid_batch) & returned_cids
        missing_cids = sorted(set(cid_batch) - returned_cids)

        all_returned_cids.update(returned_requested_cids)

        with open(
            f"outputs/batch_{batch_number:02d}.json",
            "w",
            encoding="utf-8",
        ) as file:
            json.dump(response_json, file, indent=2)

        total_batch.extend(response_json)
        

        batch_summaries.append({
            "batch": batch_number,
            "requested_count": len(cid_batch),
            "returned_requested_count": len(
                returned_requested_cids
            ),
            "missing_count": len(missing_cids),
            "missing_cids": missing_cids,
            "http_status": response.status_code,
        })

        print(
            f"Batch {batch_number:02d}: "
            f"requested={len(cid_batch)}, "
            f"returned={len(returned_requested_cids)}, "
            f"missing={len(missing_cids)}"
        )
    # Merge ATC codes into each compound object
    for compound in total_batch:
        if isinstance(compound, dict):
            cid = compound.get("cid")
            if cid is not None:
                compound["atc_code"] = cid_to_atc.get(
                    int(cid), ""
                )

    with open(
        "outputs/total_batch.json",
        "w",
        encoding="utf-8",
    ) as file:
        json.dump(total_batch, file, indent=2)

    missing_cids = sorted(
        set(requested_cids) - all_returned_cids
    )

    summary = {
        "unique_cids_in_input": len(all_cids),
        "requested_batches": NUMBER_OF_BATCHES,
        "requested_count": len(requested_cids),
        "returned_requested_count": len(all_returned_cids),
        "missing_count": len(missing_cids),
        "missing_cids": missing_cids,
        "batches": batch_summaries,
    }

    with open("outputs/summary.json", "w", encoding="utf-8") as file:
        json.dump(summary, file, indent=2)

    with open(
        "outputs/missing_cids.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as file:
        writer = csv.writer(file)
        writer.writerow(["cid", "atc_code"])

        for cid in missing_cids:
            writer.writerow([cid, cid_to_atc.get(cid, "")])

    print()
    print(f"Total requested: {summary['requested_count']}")
    print(f"Total returned:  {summary['returned_requested_count']}")
    print(f"Total missing:   {summary['missing_count']}")


if __name__ == "__main__":
    main()