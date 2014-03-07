from lxml import html
import requests

url = "http://www.rnasoft.ca/sstrand/search_results.php?select[]=Any+type&source[]=Any+source&exp_proven=Yes&select2=Greater+than+or+equal+to&first1=20&last1=&select5=No&select_duplicate=Non-redundant+sequences+only&start=0&limit=9999&sort_by=molecule%20id&order=ascending"
page = requests.get(url)
tree = html.fromstring(page.text)
tds = tree.xpath('//a[starts-with(@href, "show_results.php?molecule_ID=")]/text()')
pages = []
for id in tds:
	text = requests.get("http://www.rnasoft.ca/sstrand/show_file.php?format=Dot-parentheses&molecule_ID=" + id).text

	text = text[text.find('<textarea cols="100" rows="25" readonly align="center">'):text.find('</textarea>')]

	lines = text.split('\n')
	
	i = 6
	rna = ""
	while lines[i].strip() != "" and lines[i].strip()[0] != '(' and lines[i].strip()[0] != '.' and lines[i].strip()[0] != '[':
		rna += lines[i].strip()
		i += 1

	print(id)
	print(rna)

	sstruct = ""

	while i < len(lines):
		sstruct += lines[i].strip()
		i += 1

	print(sstruct)

	if len(rna) != len(sstruct):
		print("ERROR")
		print(str(len(rna)))
		print(str(len(sstruct)))
		print(lines)
		print(text)

		for l in lines:
			s = l.strip()
			if l not in rna and l not in sstruct:
				print(s)

		raise Exception()
