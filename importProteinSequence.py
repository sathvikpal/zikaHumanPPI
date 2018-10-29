import sys

file = sys.argv[1]

def importSequenceFromName(stringOfNames, typeOfName):
	import urllib,urllib2
	url = 'https://www.uniprot.org/uploadlists/'

	params = {'from':'ACC+ID', 'to':'ACC', 'format':'fasta', 'query':stringOfNames}

	data = urllib.urlencode(params)
	request = urllib2.Request(url, data)
	contact = "" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
	request.add_header('User-Agent', 'Python %s' % contact)
	response = urllib2.urlopen(request)
	page = response.read(2000000)
	f= open("mentha.fasta","w+")
	f.write(page)
	f.close()

f=open(file, "r")
contents =f.read()
f.close()

importSequenceFromName(contents, "")
