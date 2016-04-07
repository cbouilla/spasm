import re

# cat linbox_gauss.log | grep --invert-match "FAIL" > easy.txt

PREFIX = '/usr/share/hpac.imag.fr/Matrices/'

total_time = 0
n = 0

def parse_log(name):
	D = {}
	with open(name, 'r') as f:
		for line in f:
			match = re.search(r'(?P<filename>.+) --- (?P<time>([0-9.]+|FAIL))', line)
			assert match
			fn = match.group('filename')
			if match.group('time') == 'FAIL':
				t = -1
			else:
				t = float(match.group('time'))
			assert fn.startswith(PREFIX)
			fn = fn[len(PREFIX):]
			D[fn] = t
	return D	
		

linbox = parse_log('linbox_gauss.log')
spasm = parse_log('spasm.log')

both = { i:(linbox[i], spasm[i]) for i in linbox }

linbox_only = { m for (m, (i, j)) in both.items() if i>= 0 and j < 0}
spasm_only = { m for (m, (i, j)) in both.items() if i < 0 and j >= 0}

#print("Linbox Only : ")
for m in linbox_only:
	print(m)