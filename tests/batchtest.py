import sys

def printstuff(stuff1, stuff2):
    print(str(stuff1) + '_' + str(stuff2))

if __name__ == "__main__":
    stuff1 = sys.argv[1]
    stuff2 = sys.argv[2]

    printstuff(stuff1, stuff2)
    

