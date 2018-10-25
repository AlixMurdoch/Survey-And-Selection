import sys
import util


def importFromFile(species):
    NistFilename = util.getSpectralDataDir() + species + ' NIST.dat'

    badSymbols = [' ', '|', '[', ']', '-', '+', 'x', '?', '*', 'w', 'x']
    badSymbolsTerm = [' ', '|', '-', '+', 'x', '?', '*']

    with open(NistFilename) as NistFile:
        with open(util.getSpectralDataDir() + species + '.dat', 'w') as writeFile:
            writeFile.write('config, term, J, epsion\n')

            next(NistFile)
            row = NistFile.readline()
            
            levelIdx = row.index('L')
            termIdx = row.index('T')
            jIdx = row.index('J')

            next(NistFile)
            next(NistFile)
            oldTerm = ''
            oldConfig = ''

            for row in NistFile:
                
                config = row[0:termIdx - 3]
                term = row[termIdx-2: termIdx+5]
                j = row[jIdx-4:jIdx+5]
                level = row[levelIdx-5:levelIdx+15]
                
                for symbol in badSymbols:
                    j = j.replace(symbol, '')
                    level = level.replace(symbol, '')
                
                for symbol in badSymbolsTerm:
                    config = config.replace(symbol, '')
                    term = term.replace(symbol, '')  

                if term == '':
                    term = oldTerm
                if config == '':
                    config = oldConfig    

                
                if '<' in config:
                    openParIdx = config.index('<')
                    closeParIdx = config.index('>')
                    coreJ = config[openParIdx:closeParIdx+1]

                    config = config[:openParIdx] + str(coreJ) + config[closeParIdx+1:]
                
                if '/' in j:
                    j = float(j.split('/')[0])/float(j.split('/')[1])
                    j = str(j)
                
                #### Hack. ADD IN ENTRIES FOR LK COUPLING #######
                if '[' not in term:
                    if not row[0].isalpha(): 
                        if j.strip():
                            if level == '':
                                level = 'None'
                            writeFile.write(config + ', ')
                            writeFile.write(term + ', ')
                            writeFile.write(j + ', ')
                            writeFile.write(level + '\n')
                        else:
                            writeFile.write('\n')
                    
                oldConfig = config
                oldTerm = term
                

importFromFile(sys.argv[1])