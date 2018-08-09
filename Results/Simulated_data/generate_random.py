from random import *
from random import shuffle

bins = 63025520/3000
print bins # 21008

position = []
for i in range(500,2500):
	position.append(i*bins)

shuffle(position)


#position = []
#while len(position) < 2000:
#	num = randint(1,63025520)
#	if num not in position:
#		position.append(num)
#	else:
#		continue
#position2 = sorted(position)

#1 = deletion, 800
#2 = insertion, 800
#3 = duplication, 400

item = ["insertion"]*800
item = item + (["deletion"]*800)
item = item + (["duplication"]*400)
item_P = zip(item,position)

print isinstance(item_P[0],tuple)

list_P = []
for j in range(0,2000):
	list_P.append(list(item_P[j]))


# insertion start length
# deletion start length
# duplication start length to

for i in range(0,1600):
	size = randint(1,2001)
	list_P[i].append(size)

for v in range(1600,2000):
	size = randint(1,2001)
	list_P[v].append(size)
	to = position[v]+size  #tandam duplication
	list_P[v].append(to)

list_P = sorted(list_P,key=lambda x: x[1])

f = open("variants3.txt","a")
for thing in list_P:
	if len(thing) == 3:
		f.write(thing[0]+"	"+str(thing[1])+"	"+str(thing[2])+"\n")
	elif len(thing) == 4:
		f.write(thing[0]+"	"+str(thing[1])+"	"+str(thing[2])+"	"+str(thing[3])+"\n")
	
f.close()












































