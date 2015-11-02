
import argparse

parser = argparse.ArgumentParser(description="Compares two lists and tells how many are missing from each one and how many are found in both")
parser.add_argument("l1")
parser.add_argument("l2")
args = parser.parse_args()

itm_dict = {}
for indx, lst in enumerate([args.l1, args.l2]):
    indx += 1
    with open(lst, 'r') as IN:
        for line in IN:
            itm = line.rstrip()

            # check if item already in dict
            try:
                # if it was already found in this list, raise error
                if itm_dict[itm] >= indx:
                    raise ValueError("Duplicate name(s) in list {}. Example: '{}'".format(str(indx), itm))
                else:
                    itm_dict[itm] += indx

            # add new item to list
            except KeyError:
                itm_dict[itm] = indx

"""
parse the lists; 
in l1 but not l2 = 1
in l2 but not l1 = 2
in both = 3
"""
l1 = []
l2 = []
both = []
for k, v in itm_dict.items():
    if v == 1:
        l1.append(k)
    elif v == 2:
        l2.append(k)
    elif v == 3:
        both.append(k)
    else:
        print((k, v))
        raise Exception("This should have never happened...")

# print summaries
print("{} items in both lists:".format(str(len(both))))
for i in both:
    print(i)

print()
print("{} items in List 1 only:".format(str(len(l1))))
for i in l1:
    print(i)

print()
print("{} items in List 2 only:".format(str(len(l2))))
for i in l2:
    print(i)

