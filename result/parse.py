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
linbox = parse_log('nrow_greater_than_ncols_transpose.log')
spasm = parse_log('spasm.log')

both = { i:(linbox[i], spasm[i]) for i in linbox }

linbox_only = { m for (m, (i, j)) in both.items() if i>= 0 and j < 0}
spasm_only = { m for (m, (i, j)) in both.items() if i < 0 and j >= 0}
both_done = { m:(i,j) for (m, (i, j)) in both.items() if i>= 0 and j >= 0 }
both_fail = { m:(i,j) for (m, (i, j)) in both.items() if i < 0 and j < 0 }

print(len(linbox_only), len(spasm_only))

#linbox_faster = { m:(i, j) for (m, (i, j)) in both_done.items() if i<j}
#spasm_faster = { m:(i,j) for (m, (i, j)) in both_done.items() if i>j}

#print(len(linbox_only))
#print(len(spasm_only))
#print(len(both_done))

t_spasm = 0
t_linbox = 0
draw = 0
for m, (i,j) in both_done.items():
    #if i < j:
    #    t_linbox += 1
    #elif i > j:
    #    t_spasm += 1
    #else:
    #    draw += 1
    t_linbox += i
    t_spasm += j

print(t_linbox, t_spasm, draw)
 