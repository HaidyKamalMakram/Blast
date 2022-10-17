#Merna Hesham Mahmoud
#Haidy Kamal Makram
#Yara Amr Ahmed
#Maram Nasser
#Maram Hatem
#Yasmine Mohamed El-Gazar


from Bio.SubsMat import MatrixInfo
b62 = MatrixInfo.blosum62


#
from Bio.SubsMat import MatrixInfo
b62 = MatrixInfo.blosum62
import numpy as np

####################################################################################################################################################################################
aminoAcidSeq = 'ARNDCQEGHILKMFPSTWYV'
EnterProtein = 'PNTCGVNVDVMMERR'
print('EnterProtein= ', EnterProtein)
def Blast():

    # Open File And Read DataBase Of Protein
    OpenFile = open('D:/DataBaseOfProtein.txt', 'r')
    linesOfDB_Protein = OpenFile.readlines()

    for index, line in enumerate(linesOfDB_Protein):
        print("Protein Line {}: {}".format(index, line.strip()))

###################################################################
    ## Step [1]
    # 1st Frame
    StrFrame1 = []
    i = 0
    for i in range(len(EnterProtein)):
            StrFrame1.append(EnterProtein[i:i + 2])

    strTemp1 = []
    for indexElement in StrFrame1:
        strTemp1.append(indexElement.strip())

    np_lists_StrFrame1 = np.array(strTemp1)
    filt_StrFrame1 = []
    for i in range(0, len(strTemp1)):
        filt_i_StrFrame1 = len(np_lists_StrFrame1[i]) == 2
        filt_StrFrame1.append(filt_i_StrFrame1)

    new_list_StrFrame1 = list(np_lists_StrFrame1[filt_StrFrame1])
    strTemp1 = new_list_StrFrame1

    ####### ####### ####### ####### ####### #######
    # 2nd Frame

    StrFrame2 = []
    for i in range(len(EnterProtein)):
            StrFrame2.append(EnterProtein[i+1:i+3])

    strTemp2 = []
    for element in StrFrame2:
        strTemp2.append(element.strip())

    np_lists2_StrFrame2 = np.array(strTemp2)
    filt2_StrFrame2 = []

    for i in range(0, len(strTemp2)):
        filt_i2_StrFrame2 = len(np_lists2_StrFrame2[i]) == 2
        filt2_StrFrame2.append(filt_i2_StrFrame2)

    new_list2_StrFrame2 = list(np_lists2_StrFrame2[filt2_StrFrame2])
    strTemp2 = new_list2_StrFrame2

    ####### ####### ####### ####### ####### #######
    # 3rd Frame
    StrFrame3 = []

    for i in range(len(EnterProtein)):
            StrFrame3.append(EnterProtein[i+2:i+4])

    strTemp3 = []
    for element in StrFrame3:
        strTemp3.append(element.strip())

    np_lists3_StrFrame3 = np.array(strTemp3)
    filt3_StrFrame3 = []
    for i in range(0, len(strTemp3)):
        filt_i3_StrFrame3 = len(np_lists3_StrFrame3[i]) == 2
        filt3_StrFrame3.append(filt_i3_StrFrame3)

    new_list3_StrFrame3 = list(np_lists3_StrFrame3[filt3_StrFrame3])
    strTemp3 = new_list3_StrFrame3

    SeqFram1List = []
    for i in range(0, len(strTemp1), len(EnterProtein)-2):
        SeqFram1List.append(strTemp1[i:i + len(EnterProtein)-2])

    SeqFram2List = []
    for i in range(0, len(strTemp2), len(EnterProtein) - 2):
        SeqFram2List.append(strTemp2[i:i + len(EnterProtein) - 2])

    SeqFram3List = []
    for i in range(0, len(strTemp3), len(EnterProtein) - 2):
        SeqFram3List.append(strTemp3[i:i + len(EnterProtein) - 2])

    ####### ####### ####### ####### ####### #######
    ## Calculate Low Complexity

    lowComplexitySequence = []

    for k in range(len(SeqFram3List)):
        for i in range(len(SeqFram3List[k])):

            if i < len(SeqFram2List[k]):

                if ((b62[(SeqFram2List[k][i][0], SeqFram2List[k][i][0])]) == (b62[(SeqFram1List[k][i][0], SeqFram1List[k][i][0])]) == (b62[(SeqFram3List[k][i][0], SeqFram3List[k][i][0])])):

                    if ((b62[(SeqFram2List[k][i][1], SeqFram2List[k][i][1])]) == (b62[(SeqFram1List[k][i][1], SeqFram1List[k][i][1])]) == (b62[(SeqFram3List[k][i][1], SeqFram3List[k][i][1])])):
                        lowComplexitySequence.append(SeqFram1List[k][i])
                        lowComplexitySequence.append(SeqFram2List[k][i])
                        lowComplexitySequence.append(SeqFram3List[k][i])

        index_lowComplexitySequence = []
        for index, foundIndexInStrFrame1 in enumerate(SeqFram1List[k]):
            if foundIndexInStrFrame1 in lowComplexitySequence:
                index_lowComplexitySequence.append(index)

        for i in range(len(lowComplexitySequence)):
            SeqFram1List[k] = ['XX' if x == lowComplexitySequence[i] else x for x in SeqFram1List[k]]

#######################################################################################################################################################
    ## Step [2]
    SeqFram1ListStr = [','.join([str(elem) for elem in sublist]) for sublist in SeqFram1List]
    SeqFram1ListStrConverted = ' '.join([str(elem) for elem in SeqFram1ListStr])
    aminoAcid = list(aminoAcidSeq)
    newStrings = []
    z = 0
    for string in aminoAcid:
        newString = SeqFram1ListStrConverted.replace(",", aminoAcid[z])
        newStrings.append(newString)
        z = z + 1

    i = 0
    string1 = ''
    for i in newStrings:
        string1 = string1 + i + " "
        convertString1 = ''
        convertString1 = string1.replace(' ', '')

    i = 0
    query_wordsConvert = []
    query_seqConvert_length = int(len(convertString1))
    wordsConvert_length = query_seqConvert_length - 2
    while i < wordsConvert_length:
        query_wordsConvert.append(convertString1[i:i + 3])
        i += 3

    ConvertlistToStr = SeqFram1ListStrConverted.replace(',','')

    ConvertlistToStr2 = []
    i = 0
    for i in range(len(ConvertlistToStr)):
        ConvertlistToStr2.append(ConvertlistToStr[i:i+3])

    stringList2 = []
    i = 0
    for i in range(len(ConvertlistToStr2)):
        if i != len(ConvertlistToStr2):
            stringList2.append(ConvertlistToStr2[i][1:3])

    i = 0
    stringList2Temp = ''
    for i in stringList2:
        stringList2Temp = stringList2Temp + i + ','

    stringList2_new_strings = []
    z = 0
    for string in aminoAcid:
        s2_new_string = stringList2Temp.replace(",", "")
        stringList2_new_strings.append(s2_new_string)
        z = z + 1

    SecondString = ''
    for i in stringList2_new_strings:
        SecondString = SecondString + i + " "

    convertSecondString = ''
    convertSecondString = SecondString.replace(' ', '')

    ###############################################
    z = 0
    count = 0
    newConvertSecondString = []
    copy = stringList2Temp.replace(",", "")
    j = 0
    i = 0
    for i in range(len(aminoAcid)):
        for j in range((len(copy))):
            if len(copy[count:count + 2]):
                newConvertSecondString.insert(z, aminoAcid[i])
                newConvertSecondString.insert(z + 1, copy[count:count + 2])
                count = count + 2
                z = z + 3
        count = 0

    addAminoString2 = ''
    for i in newConvertSecondString:
        addAminoString2 = addAminoString2 + i + " "
    addAminoString2 = addAminoString2.replace(" ", "")

    i = 0
    query_wordsConvertStr2 = []
    query_seqConvert_lengthStr2 = int(len(addAminoString2))
    wordsConvert_lengthStr2 = query_seqConvert_lengthStr2 - 2
    while i < wordsConvert_lengthStr2:
        query_wordsConvertStr2.append(addAminoString2[i:i + 3])
        i += 3

    thirdString = []
    i = 0
    for i in range(len(ConvertlistToStr2)):
        if i != len(ConvertlistToStr2):
            thirdString.append(ConvertlistToStr2[i][0:3:2])

    i = 0
    ConvertthirdString = ''
    for i in thirdString:
        ConvertthirdString = ConvertthirdString + i + ','

    newConvertThirdString = []
    z = 0
    for string in aminoAcid:
        s3_new_string = ConvertthirdString.replace(",", "")
        newConvertThirdString.append(s3_new_string)
        z = z + 1

    string3 = ''
    i = 0
    for i in newConvertThirdString:
        string3 = string3 + i + " "

    convertString3 = ''
    convertString3 = string3.replace(' ', '')

    zzz = 1
    count3 = 0
    newConvertString3 = []
    copy3 = ConvertthirdString.replace(",", "")
    i = 0
    j = 0
    for i in range(len(aminoAcid)):
        for j in range((len(copy3))):
            if len(copy3[count3:count3 + 2]):
                newConvertString3.insert(zzz, copy3[count3:count3 + 1])
                newConvertString3.insert(zzz + 1, aminoAcid[i])
                newConvertString3.insert(zzz + 2, copy3[count3 + 1:count3 + 2])
                count3 = count3 + 2
                zzz = zzz + 3
        count3 = 0

    stw3 = ''
    for i in newConvertString3:
        stw3 = stw3 + i + " "
    stw3 = stw3.replace(" ", "")

    ###############################################

    i = 0
    query_wordsConvert3 = []
    query_seqConvert_lengthStr3 = int(len(stw3))
    wordsConvert_lengthStr3 = query_seqConvert_lengthStr3 - 2
    while i < wordsConvert_lengthStr3:
        query_wordsConvert3.append(stw3[i:i + 3])
        i += 3

#########################################################################################################################
    #### for Frame 1
    i = 0
    j = 0
    TotalSum = 0
    for i in range(len(query_wordsConvert)):
        for j in range(len(query_wordsConvert[0])):
            if j < len(query_wordsConvert[0]):
                TotalSum = TotalSum + b62[(query_wordsConvert[i][j], query_wordsConvert[i][j])]
        if TotalSum < 13:
            del (i)
        TotalSum = 0

    ### fro Frame 2
    i = 0
    j = 0
    TotalSum2 = 0
    for i in range(len(query_wordsConvertStr2)):
        for j in range(len(query_wordsConvertStr2[0])):
            if j < len(query_wordsConvertStr2[0]):
                TotalSum2 = TotalSum2 + b62[(query_wordsConvertStr2[i][j], query_wordsConvertStr2[i][j])]
        if TotalSum2 < 13:
            del (i)
        TotalSum2 = 0

    ### for Line 3
    i = 0
    j = 0
    TotalSum3 = 0
    for i in range(len(query_wordsConvert3)):
        for j in range(len(query_wordsConvert3[0])):
            if j < len(query_wordsConvert3[0]):
                TotalSum3 = TotalSum3 + b62[(query_wordsConvert3[i][j], query_wordsConvert3[i][j])]
        if TotalSum3 < 13:
            del(i)
        TotalSum3 = 0

    i = 0
    for i in range(len(EnterProtein)):
        if i < len(EnterProtein):
            HighScoreWords1 = [item for item in query_wordsConvert if(item in EnterProtein)]
    HighScoreWords1 = list(dict.fromkeys(HighScoreWords1))
    print('HighScoreWords1: ', HighScoreWords1)

    i = 0
    for i in range(len(EnterProtein)):
        if i < len(EnterProtein):
            HighScoreWords2 = [item for item in query_wordsConvertStr2 if (item in EnterProtein)]
    HighScoreWords2 = list(dict.fromkeys(HighScoreWords2))
    print('HighScoreWords2: ', HighScoreWords2)

    i = 0
    for i in range(len(EnterProtein)):
        if i < len(EnterProtein):
            HighScoreWords3 = [item for item in query_wordsConvert3 if (item in EnterProtein)]

    #remove repeat
    HighScoreWords3 = list(dict.fromkeys(HighScoreWords3))
    print('HighScoreWords3: ', HighScoreWords3)

    StopScore = 3
    totalScoreForEachHigh = 0
    i = 0
    j = 0

    Match1 = []
    i = 0
    j = 0
    k =0
    for i in range(len(linesOfDB_Protein)):
        if i < len(linesOfDB_Protein):
            for j in range(len(linesOfDB_Protein[i])):
                if j < len(linesOfDB_Protein):
                    for k in range(len(HighScoreWords1)):
                        if linesOfDB_Protein[i][j : j + 3] == HighScoreWords1[k]:
                            if j == 0:
                                totalScoreForEachHigh = b62[(linesOfDB_Protein[i][j : j + 3], linesOfDB_Protein[i][j : j + 3])] + b62(HighScoreWords1[k])
                                Match1.append(HighScoreWords1[k]+linesOfDB_Protein[j+4])
                                #if totalScoreForEachHigh + linesOfDB_Protein[i][]
    print('Match1: ', Match1)

    Match2 = []
    i = 0
    j = 0
    for i in range(len(linesOfDB_Protein)):
        if i < len(linesOfDB_Protein):
            for j in range(len(linesOfDB_Protein[i])):
                if j < len(linesOfDB_Protein):
                    Match2 = [ele for ele in HighScoreWords2 if (ele in linesOfDB_Protein[i][j:j + 3])]

    print('Match1: ', Match2)

    Match3 = []
    i = 0
    j = 0
    for i in range(len(linesOfDB_Protein)):
        if i < len(linesOfDB_Protein):
            for j in range(len(linesOfDB_Protein[i])):
                if j < len(linesOfDB_Protein):
                    Match3 = [ele for ele in HighScoreWords3 if (ele in linesOfDB_Protein[i][j:j + 3])]

    print('Match3: ', Match3)




Blast()