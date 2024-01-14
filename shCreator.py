filePrefix = r'/home1/flippenc/shFiles'
for i in range(1,21):
    for j in [1000,2000,5000]:
        for k in [200,500]:
            with open(f'{filePrefix}/paramNorm-{i/10}-{j}-{k}.sh','w') as shFile:
                shFile.write('#!/bin/bash\n')
                shFile.write('module load matlab/R2021a\n\n')
                shFile.write(f'matlab -nodisplay -nodesktop -r "parameterizedGenModelElliptical 0 0 0 1 0.01 {i/10} 1 0.5 {j} {k}"')
for i in range(3,11):
    for j in [1000,2000,5000]:
        for k in [200,500]:
            with open(f'{filePrefix}/paramNorm-{i}-{j}-{k}-{l}.sh','w') as shFile:
                shFile.write('#!/bin/bash\n')
                shFile.write('module load matlab/R2021a\n\n')
                shFile.write(f'matlab -nodisplay -nodesktop -r "parameterizedGenModelElliptical 0 0 0 1 0.01 {i} 1 0.5 {j} {k}"')
