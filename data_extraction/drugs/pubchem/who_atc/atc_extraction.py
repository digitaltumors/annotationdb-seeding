import requests
import json
import csv

atc_data = {
	"Annotations": []
}
for i in range(6):
	response = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/JSON/?source=WHO+Anatomical+Therapeutic+Chemical+%28ATC%29+Classification&heading=ATC+Code&page={i+1}&heading_type=Compound&response_type=display")
	response_json = response.json()
	atc_data["Annotations"].extend(response_json.get('Annotations').get('Annotation'))
	# 	atc_data["Annotations"] = atc_data["Annotations"] + response_json.get('Annotations').get('Annotation')

with open("atc_data.json", "w") as f:
	json.dump(atc_data, f, indent=4)

atc_codes = []
existing_codes = set()

for annotation in atc_data["Annotations"]:
	for atc_entry in annotation["Data"][0]["Value"]["StringWithMarkup"]:
		code, description = atc_entry["String"].split(" - ", 1)

		if code not in existing_codes:
			atc_codes.append({
				"code": code,
				"description": description
			})

			existing_codes.add(code)

atc_codes.sort(key=lambda x: x["code"]) # sort to easily read code csv


with open("atc_codes.csv", "w", newline="") as f:
	writer = csv.DictWriter(
		f,
		fieldnames=["code", "description"]
	)

	writer.writeheader()
	writer.writerows(atc_codes)

	
	