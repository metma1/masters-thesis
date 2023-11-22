import parameters as p


def ConcentrationInit():
    OxyConcList = ''
    NutriConcList = ''

    for x in range(p.Xdim+1):
        for y in range(p.Ydim+1):
            for z in range(p.Zdim+1):
                if y < p.Agar_Thickness:
                    OxyConcList += str(x) + '\t' + str(y) + '\t'  + str(z) + '\t' + '0' + '\n'
                    NutriConcList += str(x) + '\t' + str(y) + '\t'  + str(z) + '\t' + str(p.Initial_Nutrient) + '\n'
                else:
                    OxyConcList += str(x) + '\t' + str(y) + '\t'  + str(z) + '\t' + str(p.Initial_Oxygen) + '\n'
                    NutriConcList += str(x) + '\t' + str(y) + '\t'  + str(z) + '\t' + '0' + '\n'


    with open("InitOxygenContent.txt", 'w') as OxyF:
        OxyF.write(OxyConcList)
    OxyF.close()

    with open("InitNutrientContent.txt", 'w') as NutriF:
        NutriF.write(NutriConcList)
    NutriF.close()

ConcentrationInit()