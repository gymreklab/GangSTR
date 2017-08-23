read_dict_pyt = {}
read_dict_cpp = {}
with open ('python_span.txt', 'r') as pyt:
	with open ('cpp_span.txt', 'r') as cpp:
		for record in pyt:
			spl = record.split('\t')
			if spl[0] != '' and len(spl) == 3:
				if spl[0] not in read_dict_pyt:
					read_dict_pyt[spl[0]] = spl[1].strip()
					# read_dict_pyt[spl[0]] = spl[1] + spl[2].strip()
				# else:
				# 	read_dict[spl[0]] = read_dict[spl[0]] + " " + spl[1] + spl[2].strip()
		for record in cpp:
			spl = record.split('\t')
			if spl[0] != '' and len(spl) == 2:
				if spl[0] not in read_dict_cpp:
					read_dict_cpp[spl[0]] = spl[1].strip()
					# read_dict_cpp[spl[0]] = spl[1] + spl[2].strip()

		for key in read_dict_pyt:
			if key in read_dict_cpp:
				if read_dict_pyt[key] != read_dict_cpp[key]:
					pass
					# print key, '\t', read_dict_pyt[key], '\t -> CPP: ', read_dict_cpp[key]
			else:
				print key, '\tNot in CPP'
		for key in read_dict_cpp:
			if key not in read_dict_pyt:
				print key, '\tNot in Python'	
