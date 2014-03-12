#Process US Census data to give a yearly population estimate for 2000-2011 
# which matches the list of cities used to produce the distance matrix. 

census2009 = open("../scriptData/census/SUB-EST2009-04-19.csv")

def getCitiesObject():
    '''Read the list_of_cities_in_iowa.txt file into memory'''
    outList = []
    citiesFile = open("../scriptData/OSM/list_of_cities_in_iowa.txt", "r")
    for i in citiesFile.readlines():
        outList.append(i.strip("\r\n"))
    citiesFile.close()
    return(outList) 

def readCensus2011():
    ''' Read the most recent census data into memory, county information is available here if that would be useful in the future'''
    outList = []
    census2011 = open("../scriptData/census/SUB-EST2011_19.csv", "r")
    firstLine = True
    for line in census2011.readlines():
        if firstLine:
            firstLine = False
            continue
        lineSplit = line.split(",")
        cityName = lineSplit[6]
        census2010Pop = lineSplit[8]
        census2011Pop = lineSplit[9]
        if  ("city" in cityName) and (int(lineSplit[0]) == 162):
            cityName = cityName.split(" ")

            outList.append((" ".join(cityName[0:(len(cityName)-1)]), census2010Pop, census2011Pop))
    return(outList )

def readCensus2009():
    ''' Read census data for 2000 to 2009 into memory, county information is available here if that would be useful in the future'''
    outList = []
    census2009 = open("../scriptData/census/SUB-EST2009-04-19_v2.csv", "r")
    lineIdx = -1
    for line in census2009.readlines():
        lineIdx += 1
        if lineIdx < 1:
            # Skip the header
            continue

        line = line.strip("\r\n")
        lineSplit = line.split(",")
        cityName = lineSplit[0].replace(" city", "")[1:]
        jul09 = lineSplit[1]
        jul08 = lineSplit[2]
        jul07 = lineSplit[3]
        jul06 = lineSplit[4]
        jul05 = lineSplit[5]
        jul04 = lineSplit[6]
        jul03 = lineSplit[7]
        jul02 = lineSplit[8]
        jul01 = lineSplit[9]
        jul00 = lineSplit[10]
        outList.append((cityName, jul00, jul01, jul02, jul03, jul04, jul05, jul06, jul07, jul08, jul09))
    return(outList)

listOfCities = getCitiesObject()
census2009 = readCensus2009()
census2011 = readCensus2011()

def which(listVar):
    idx = 0
    for item in listVar:
        if item:
            return(idx)
        idx += 1

def combineCensusData():
    outList = []
    noFind11 = []
    noFind09 = []
    error = []
    for city in listOfCities:
        idx09 = [x == city for x in [rw[0] for rw in census2009]]
        idx11 = [x == city for x in [rw[0] for rw in census2011]]
        evt = 0
        if (sum(idx09) == 0):
            evt += 1
            noFind09.append(city)
        if (sum(idx11) == 0):
            evt += 1
            noFind11.append(city)
        if ((sum(idx11) == 1) and (sum(idx09) == 1)):
            evt += 1
            c09 = census2009[which(idx09)]
            c11 = census2011[which(idx11)]
            outList.append(([city] + [x for x in c09[1:]] + [x for x in c11[1:]]))
        if evt == 0:
            error.append(city)
    print("Success: " + str(len(outList)))
    print("Missing 09: " + str(len(noFind09)))
    print("Missing 11: " + str(len(noFind11)))
    print("Error: " + str(len(error)))


    return({"outList":outList, "noFind11":noFind11, "noFind09":noFind09, "error":error})

MatchedData = combineCensusData()

# The following places were unable to be matched
# ['Ira', 'Kent', 'Littleport', 'Placid', 'Vedic City']

outFile = open("../scriptData/census/ProcessedPopulationData.csv", "w")
outFile.write("City;y00;y01;y02;y03;y04;y05;y06;y07;y08;y09;y10;y11\n")
for i in MatchedData["outList"]:
    outFile.write(";".join([str(x) for x in i]) + "\n")
outFile.close()





