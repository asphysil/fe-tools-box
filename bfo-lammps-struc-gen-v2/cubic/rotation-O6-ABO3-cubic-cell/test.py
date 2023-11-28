# importing the module
from collections import Counter
 
# making a list
list = [1, 1, 2, 3, 4, 5,
        6, 7, 9, 2, 3, 4, 8]
 
# instantiating a Counter object
ob = Counter(list)
 
# Counter.items()
items = ob.items()
 
print("The datatype is "
      + str(type(items)))
 
# displaying the dict_items
print(items)
 
# iterating over the dict_items
#for i in items:
#    print(i)