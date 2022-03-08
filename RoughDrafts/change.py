import io
import sys
import pandas as pd

def main():

    #save the original standard output for outputting to terminal
    stdout_fileno = sys.stdout

    #redirect standard out to our file
    sys.stdout = open('./data_set.csv', 'w')

    data = open('../GSE64881_segmentation_at_30000bp.passqc.multibam.txt', 'r')
    for line in data:
        newline = ', '.join(line.split('\t'))
        sys.stdout.write(newline)

    #revert stdout to original
    sys.stdout.close()
    sys.stdout = stdout_fileno

    print('converted')

    return

main()
