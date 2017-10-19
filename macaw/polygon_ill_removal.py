#Polygonal Illegal Path Remover.

#This method will take a list, of some numbers, and "+" and "-". If
#this given path is not the most efficient, it is considered an illegal path.

#This method will get rid of the illegal parts of the path, and make
#the path legal.


def type(s):
        if s == "+" or s == "-":
            return 1
        elif s == "++" or s == "--":
            return 2
        elif s == "C" or s == "A":
            return 3
        else:
            return 0

def Polygon_Ill_Remover (path_list):
    list1 = path_list;
    i = 0

    while(i < len(list1) - 2):
        print(list1);
        if (list1[i] == "+" or list1[i+2] == "-") and  (list1[i] == list1[i+2]) and (type(list1[i+1]) == 0):
            del list1[i];
            del list1[i];
            del list1[i];
            i -= 2
        i += 1

# This is an example: Polygon_Ill_Remover([1,"+","+",3,"-","-",4,"-","+",3,"+","-",3,"-","-",1]);
