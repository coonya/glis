from ConfigParser import SafeConfigParser
import sys, optparse

def read_config(config_file):
	params = {}
	parser = SafeConfigParser()
	parser.read(config_file)

	for name,value in parser.items('common'):
		params[name] = value
	

	return params


def check_dir(dir):
	if os.path.exists(dir):
		pass
	else:
		try:
			os.mkdir(dir)
		except:
			pass



def panel(num):
	if num == '1':
		panel = 'OncoPanel_v2'
	elif num == '2':
		panel = 'AMCv1_C'
	elif num == '3':
		panel = 'AMCv1_S'
	elif num == '4':
		panel = 'AMCv2'
	elif num == '5':
		panel = 'Tier2'
	elif num == '6':
		panel = 'AMCv3C'
	elif num == '7':
		panel = 'AMCv3S'
	elif num == '8':
		panel = 'HEMEv1'
	
	return panel
