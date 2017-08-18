read_dict_unk = {}
read_dict_rch = {}
with open ('output_just_unknwn_recheck.txt', 'r') as unk:
	with open ('output_with_recheck.txt', 'r') as rch:
		for record in unk:
			spl = record.split('\t')
			if spl[0] != '' and len(spl) == 4:
				if spl[0] not in read_dict_unk:
					# read_dict_unk[spl[0]] = spl[1]
					read_dict_unk[spl[0]] = spl[1] + spl[2].strip()
				# else:
				# 	read_dict[spl[0]] = read_dict[spl[0]] + " " + spl[1] + spl[2].strip()
		for record in rch:
			spl = record.split('\t')
			if spl[0] != '' and len(spl) == 4:
				if spl[0] not in read_dict_rch:
					read_dict_rch[spl[0]] = spl[1] + spl[2].strip()

		for key in read_dict_rch:
			if key in read_dict_unk:
				if read_dict_unk[key] != read_dict_rch[key]:
					print key, '\t', read_dict_unk[key], '\t -> recheck_all: ', read_dict_rch[key]
			else:
				print key, ' Not present when only checking UNKNOWN'
		for key in read_dict_unk:
			if key not in read_dict_rch:
				print key, ' Not present when rechecking all.'	
