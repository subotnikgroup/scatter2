
# Specialized script for preprocesing the compilation


import sys
src = sys.argv[1:-1]
dst = sys.argv[-1]

# collect class names
class_names = []
for s in src:
    class_names += [line.strip().split()[1].split(':')[0] for line in open(s, 'r') 
        if line.strip().startswith('class ')]
for line in open(dst, 'r'):
    if line.strip().startswith('preprocfac'):
        cname = line.split()[1].split(';')[0]
        for c in class_names:
            print('    %s (strcmp(%s, %s::get_name()) == 0){return new %s();}' % 
            ('if' if c == class_names[0] else 'else if', cname, c, c))
        print('    else {throw std::runtime_error("Potential not found\\n");}')
    else:
        print(line.strip('\n'))

