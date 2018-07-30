# Quick script to compare analyses

# open up the files
meAllFile = open("AllCutsCutEvents.txt", 'r')
stefaniaAllFile = open("StefaniaAll.txt", 'r')

# read into lists
meAll = []
stefaniaAll = []
for line in meAllFile:
    meAll.append(int(line.strip()))
for line in stefaniaAllFile:
    stefaniaAll.append(int(line.strip())-1)

meNotStefaniaAll = []
stefaniaNotMeAll = []

for me in meAll:
    if me not in stefaniaAll:
        meNotStefaniaAll.append(me)
for stefania in stefaniaAll:
    if stefania not in meAll:
        stefaniaNotMeAll.append(stefania)

print "All cuts:"
print "Me not Stefania:"
print meNotStefaniaAll
print "Stefania not me:"
print stefaniaNotMeAll
