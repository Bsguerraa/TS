import string

class Classe:

    def __init__(self, dominios, sequencias, output):
        self.dominios = dominios
        self.sequencias = sequencias
        self.output = output
        self.domnomes = ['HAMP', 'HisKA', 'HATPase_c', 'Response_reg', 'Sigma54_activat', 'HTH_8']

    def range(self, nomeseq, nomedom):
        arqdom = open(self.dominios, 'r').readlines()
        contador = 0
        primeiro = ""
        segundo = ""
        intervalo = []

        for line in arqdom:
            if nomeseq in line:
                if nomedom in line:
                    charinicial = line.index(" ") + 1
                    while contador < 6:
                        if line[charinicial] == " ":
                            if contador != 0:
                                break
                            else:
                                charinicial = charinicial + 1
                                contador = 0
                        else:
                            contador = contador + 1
                            primeiro = primeiro + str(line[charinicial])
                            charinicial = charinicial + 1

                    charinicial = charinicial + 1
                    contador = 0

                    while contador < 4:
                        if line[charinicial] == " ":
                            if contador != 0:
                                break
                            else:
                                charinicial = charinicial + 1
                                contador = 0
                        else:
                            contador = contador + 1
                            segundo = segundo + str(line[charinicial])
                            charinicial = charinicial + 1
        intervalo.append(primeiro)
        intervalo.append(segundo)
        return intervalo

    def nomes(self):
        contador = 0
        name = ""
        listanomes = []
        arq = open(self.sequencias, 'r').read()
        namestart = False
        while contador < arq.__len__():

            if namestart:
                name = name + arq[contador]
                if arq[contador] == " ":
                    namestart = False
                    listanomes.append(name)
                    name = ""

            if arq[contador] == ">":
                namestart = True

            contador += 1

        print(arq.count(">"), "sequências encontradas")
        return listanomes

    def listasequencias(self):
        import re
        onlySequences = []
        arqseq = open(self.sequencias, 'r').read()
        regex = re.compile(r'^>([^\n\r]+)[\n\r]([A-Z\n\r]+)', re.MULTILINE)
        matches = [m.groups() for m in regex.finditer(arqseq)]

        for m in matches:
            onlySequences.append(m[1])
        return onlySequences

    def transformar_lowercase(self, value, lower, upper):
        return value[:lower-1] + value[lower:upper].lower() + value[upper+1:]

    def retornarlinha(self):
        arq = open(self.sequencias, 'r').readlines()
        linhas_sem_sequencias = []
        for line in arq:
            if ">" in line:
                linhas_sem_sequencias.append(line)
        return linhas_sem_sequencias

    def highlight(self):
        arq = open(self.sequencias, 'r').read()
        novoarquivo = open(self.output, 'w')
        a = 0

        linhas_sem_sequencias = self.retornarlinha()
        listasequencias = self.listasequencias()
        listanomes = self.nomes()

        while a < arq.count(">"):

            HAMP = False
            HisKA = False
            HATPase_c = False
            Response_reg = False
            Sigma54_activat = False
            HTH_8 = False

            novoarquivo.write(linhas_sem_sequencias[a])

            while not(Response_reg):
                x = self.range(listanomes[a], "Response_reg")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                Response_reg = True

            while not(Sigma54_activat):
                x = self.range(listanomes[a], "Sigma54_activat")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                Sigma54_activat = True

            while not(HTH_8):
                x = self.range(listanomes[a], "HTH_8")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HTH_8 = True

            while not(HAMP):
                x = self.range(listanomes[a], "HAMP")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HAMP = True

            while not(HisKA):
                x = self.range(listanomes[a], "HisKA")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HisKA = True

            while not (HATPase_c):
                x = self.range(listanomes[a], "HATPase_c")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HATPase_c = True

            novoarquivo.write(listasequencias[a])

            a+=1

    def rfinal(self):
        arq = open(self.sequencias, 'r').read()
        novoarquivo = open(self.output, 'w')
        arqdom1 = open('response_reg.fasta', 'w')
        arqdom2 = open('sigma54.fasta', 'w')
        arqdom3 = open('hth.fasta', 'w')
        a = 0
        b = 0
        linhas_sem_sequencias = self.retornarlinha()
        listasequencias = self.listasequencias()
        listanomes = self.nomes()

        while b < arq.count(">"):
            x = self.range(listanomes[b], "Response_reg")
            if x[0] != '':
                linhas_sem_sequencias[b] = linhas_sem_sequencias[b][:(linhas_sem_sequencias[b].rfind('Response_reg')) + 12] + ' ' + x[0] + ' a ' + x[1]  + linhas_sem_sequencias[b][(linhas_sem_sequencias[b].rfind('Response_reg')) + 12:]
            x = self.range(listanomes[b], "activat")
            if x[0] != '':
                linhas_sem_sequencias[b] = linhas_sem_sequencias[b][:(linhas_sem_sequencias[b].rfind('activat')) + 7] + ' ' + x[
                    0] + ' a ' + x[1] + linhas_sem_sequencias[b][(linhas_sem_sequencias[b].rfind('activat')) + 7:]
            x = self.range(listanomes[b], "HTH_8")
            if x[0] != '':
                linhas_sem_sequencias[b] = linhas_sem_sequencias[b][:(linhas_sem_sequencias[b].rfind('HTH_8')) + 5] + ' ' + x[
                    0] + ' a ' + x[1] + linhas_sem_sequencias[b][(linhas_sem_sequencias[b].rfind('HTH_8')) + 5:]
            x = self.range(listanomes[b], "HAMP")
            if x[0] != '':
                linhas_sem_sequencias[b] = linhas_sem_sequencias[b][:(linhas_sem_sequencias[b].rfind('HAMP')) + 4] + ' ' + x[
                    0] + ' a ' + x[1] + linhas_sem_sequencias[b][(linhas_sem_sequencias[b].rfind('HAMP')) + 4:]
            x = self.range(listanomes[b], "HisKA")
            if x[0] != '':
                linhas_sem_sequencias[b] = linhas_sem_sequencias[b][:(linhas_sem_sequencias[b].rfind('HisKA')) + 5] + ' ' + x[
                    0] + ' a ' + x[1] + linhas_sem_sequencias[b][(linhas_sem_sequencias[b].rfind('HisKA')) + 5:]
            x = self.range(listanomes[b], "HATPase_c")
            if x[0] != '':
                linhas_sem_sequencias[b] = linhas_sem_sequencias[b][:(linhas_sem_sequencias[b].rfind('HATPase_c')) + 9] + ' ' + x[
                    0] + ' a ' + x[1] + linhas_sem_sequencias[b][(linhas_sem_sequencias[b].rfind('HATPase_c')) + 9:]

            b += 1

        while a < arq.count(">"):

            HAMP = False
            HisKA = False
            HATPase_c = False
            Response_reg = False
            Sigma54_activat = False
            HTH_8 = False

            novoarquivo.write(linhas_sem_sequencias[a])
            arqdom1.write(linhas_sem_sequencias[a])
            arqdom2.write(linhas_sem_sequencias[a])
            arqdom3.write(linhas_sem_sequencias[a])

            transf = list(listasequencias)

            while not (Response_reg):
                x = self.range(listanomes[a], "Response_reg")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                regg = self.transformar_lowercase(transf[a], int(x[0]), int(x[1]))

                string1 = regg
                list1 = list(string1)
                new_list = []
                for i in list1:
                    if i.isupper():
                        i = ""
                        new_list.append(i)
                    else:
                        new_list.append(i)

                regg = ''.join(new_list)
                arqdom1.write(regg)

                Response_reg = True

            while not (Sigma54_activat):
                x = self.range(listanomes[a], "Sigma54_activat")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                sigm = self.transformar_lowercase(transf[a], int(x[0]), int(x[1]))

                string1 = sigm
                list1 = list(string1)
                new_list = []
                for i in list1:
                    if i.isupper():
                        i = ""
                        new_list.append(i)
                    else:
                        new_list.append(i)

                sigm = ''.join(new_list)
                arqdom2.write(sigm)

                Sigma54_activat = True

            while not (HTH_8):
                x = self.range(listanomes[a], "HTH_8")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                hthh = self.transformar_lowercase(transf[a], int(x[0]), int(x[1]))
                string1 = hthh
                list1 = list(string1)
                new_list = []
                for i in list1:
                    if i.isupper():
                        i = ""
                        new_list.append(i)
                    else:
                        new_list.append(i)

                hthh = ''.join(new_list)
                arqdom3.write(hthh)
                HTH_8 = True

            while not (HAMP):
                x = self.range(listanomes[a], "HAMP")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HAMP = True

            while not (HisKA):
                x = self.range(listanomes[a], "HisKA")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HisKA = True

            while not (HATPase_c):
                x = self.range(listanomes[a], "HATPase_c")
                if x[0] == "":
                    break
                listasequencias[a] = self.transformar_lowercase(listasequencias[a], int(x[0]), int(x[1]))
                HATPase_c = True

            string1 = listasequencias[a]
            list1 = list(string1)
            new_list = []
            for i in list1:
                if i.isupper():
                    i = ""
                    new_list.append(i)
                else:
                    new_list.append(i)

            listasequencias[a] = ''.join(new_list)
            novoarquivo.write(listasequencias[a])
            a += 1
