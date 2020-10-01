#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <time.h>

#define NUMALPHAS 6

using namespace std;

typedef struct {
    int id;
    int numSalas;
    int largura;
    int *fluxo = nullptr;
} Sala;

vector<Sala>* carregaInstancia(string nomeArquivo) {
    vector<Sala> *salas = new vector<Sala>;
    int numSalas;
    int *fluxo;
    ifstream arquivoEntrada;
    string str;
    stringstream ss;
    arquivoEntrada.open(nomeArquivo);

    if(arquivoEntrada.is_open()) {
        getline(arquivoEntrada, str);
        ss << str;
        ss >> numSalas;

        getline(arquivoEntrada, str);
        ss.clear();
        ss << str;
        for(int i=0; i<numSalas; i++) {
            getline(ss, str, ',');
            Sala aux;
            aux.id = i;
            aux.numSalas = numSalas;
            aux.largura = stoi(str);
            salas->push_back(aux);
        }

        for(int i=0; i<numSalas; i++) {
            getline(arquivoEntrada, str);
            ss.clear();
            ss << str;
            fluxo = new int[numSalas];
            for(int j=0; j<numSalas; j++) {
                getline(ss, str, ',');
                fluxo[j] = stoi(str);
            }
            salas->at(i).fluxo = fluxo;
        }
    } else
        cout << "Arquivo nao Encontrado" << endl;
    return salas;
}

void heapfy(vector<double> *custosCandidatos, vector<int> *idCandidatos, int i, int n) {
    int maior = i;
    int esquerda= 2*i +1;
    int direita = 2*i +2;
    if(esquerda < n && custosCandidatos->at(esquerda) < custosCandidatos->at(maior))
        maior = esquerda;
    if(direita < n && custosCandidatos->at(direita) < custosCandidatos->at(maior))
        maior = direita;
    if(maior != i) {
        swap(custosCandidatos->at(i), custosCandidatos->at(maior));
        swap(idCandidatos->at(i), idCandidatos->at(maior));
        heapfy(custosCandidatos, idCandidatos, maior, n);
    }
}

void heapSort(vector<double> *custosCandidatos, vector<int> *idCandidatos) {
    double aux;
    int n = custosCandidatos->size();
    for(int i=custosCandidatos->size()/2 - 1; i>=0; i--)
        heapfy(custosCandidatos, idCandidatos, i, n);
    for(int i=custosCandidatos->size()-1; i>=0; i--) {
        swap(custosCandidatos->at(i), custosCandidatos->at(0));
        swap(idCandidatos->at(i), idCandidatos->at(0));
        heapfy(custosCandidatos, idCandidatos, 0, i);
    }
}

double calculaCusto(vector<Sala> *salas, double *corredor) {
    double custo = 0;
    for(int i=0; i<salas->size(); i++) {
        for(int j=i+1; j<salas->size(); j++)
                custo += fabs(fabs(corredor[i]) - fabs(corredor[j]))*salas->at(i).fluxo[j];
    }
    return custo;
}

double estimaCusto(vector<Sala> *salas, double custo, double *corredor, double posCandidato, int idCandidato) {
    for(int i=0; i<salas->at(0).numSalas; i++) {
        if(corredor[i] != 0)
            custo += fabs(fabs(corredor[i]) - fabs(posCandidato))*salas->at(i).fluxo[idCandidato];
    }
    return custo;
}

void apresentaSolucao(vector<Sala> *salas, double* corredor) {
    cout << "Custo: " << calculaCusto(salas, corredor) << endl;
}

int selecionaOrigem(vector<Sala> *salas) {
    int numSalas = salas->size();
    int idMenorFluxo;
    double somatorioFluxos;
    double menorFluxo = INFINITY;

    for(int i=0; i<numSalas; i++) {
        somatorioFluxos = 0;
        for(int j=0; j<numSalas; j++)
            somatorioFluxos += salas->at(i).fluxo[j];
        if(somatorioFluxos < menorFluxo) {
            menorFluxo = somatorioFluxos;
            idMenorFluxo = i;
        }
    }
    return idMenorFluxo;
}

double* buscaLocalSwap(vector<Sala> *salas, double* corredor) {
    int numSalas = salas->size();
    double* novoCorredor = new double[numSalas];

    for(int i=0; i<numSalas; i++)
        novoCorredor[i] = corredor[i];

    for(int i=0; i<numSalas; i++) {
        for(int j=1; j<numSalas; j++) {

        }
    }
}

double* guloso(vector<Sala> *salas) {
    int numSalas = salas->size();
    int numSalasInseridas = 0;
    int idAtual;
    int idMaiorFluxo;
    int idUltimaAdd;
    double maiorFluxo;
    double *corredor = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    double larguraMaiorFluxo;
    double custo = 0;
    double custoAux;
    bool *salasInseridas = new bool[numSalas];

    for(int i=0; i<numSalas; i++) {
        corredor[i] = 0;
        salasInseridas[i] = false;
    }

    idAtual = selecionaOrigem(salas);
    espOcupadoSup = salas->at(idAtual).largura;
    corredor[idAtual] = espOcupadoSup/2;
    salasInseridas[idAtual] = true;
    numSalasInseridas++;

    idUltimaAdd = idAtual;
    while(numSalasInseridas < numSalas) {
        maiorFluxo = 0;
        for(int i=0; i<numSalas; i++) {
            if(!salasInseridas[i]) {
                custoAux = fabs((fabs(corredor[idUltimaAdd]) - espOcupadoMenor + salas->at(i).largura/2)) * salas->at(idUltimaAdd).fluxo[i];
                
                if(custoAux > maiorFluxo) {
                    idMaiorFluxo = i;
                    maiorFluxo = custoAux;
                    larguraMaiorFluxo = salas->at(i).largura;
                }
            }
        }
        
        if(espOcupadoSup < espOcupadoInf) {
            corredor[idMaiorFluxo] = espOcupadoSup + larguraMaiorFluxo/2;
            espOcupadoSup += larguraMaiorFluxo;
        } else {
            corredor[idMaiorFluxo] = (espOcupadoInf + larguraMaiorFluxo/2)*-1;
            espOcupadoInf += larguraMaiorFluxo;
        }
        if(espOcupadoSup < espOcupadoInf)
            espOcupadoMenor = espOcupadoSup;
        else
            espOcupadoMenor = espOcupadoInf;
        
        salasInseridas[idMaiorFluxo] = true;
        numSalasInseridas++;

        idUltimaAdd = idMaiorFluxo;
    }

    delete[] salasInseridas;
    return corredor;
}

double* auxGulosoRandomizado(int seed, float alpha, vector<Sala> *salas) {
    srand(seed);
    int numSalas = salas->size();
    int numSalasInseridas = 0;
    int aleatorio;
    int idSala;
    int idUltimaAdd;
    double *corredor = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    vector<int> *idCandidatos = new vector<int>;
    vector<double> *fluxoCandidatos = new vector<double>;
    bool *salasInseridas = new bool[numSalas];

    for(int i=0; i<numSalas; i++) {
        corredor[i] = 0;
        salasInseridas[i] = false;
        idCandidatos->push_back(i);
        fluxoCandidatos->push_back(0);
    }

    idSala = selecionaOrigem(salas);
    espOcupadoSup = salas->at(idSala).largura;
    corredor[idSala] = espOcupadoSup/2;
    salasInseridas[idSala] = true;
    numSalasInseridas++;

    idUltimaAdd = idSala;
    while(numSalasInseridas < numSalas) {
        for(int i=0; i<numSalas; i++) {
            idCandidatos->at(i) = i;
            if(!salasInseridas[i])
                fluxoCandidatos->at(i) = fabs((fabs(corredor[idUltimaAdd]) - espOcupadoMenor + salas->at(i).largura/2)) * salas->at(idUltimaAdd).fluxo[i];
            else
                fluxoCandidatos->at(i) = 0;
        }
        heapSort(fluxoCandidatos, idCandidatos);
        
        aleatorio = rand()%(int)ceil(alpha*(numSalas-numSalasInseridas));
        fluxoCandidatos->at(aleatorio);
        idSala = idCandidatos->at(aleatorio);
        if(espOcupadoSup < espOcupadoInf) {
            corredor[idSala] = espOcupadoSup + salas->at(idSala).largura/2;
            espOcupadoSup += salas->at(idSala).largura;
        } else {
            corredor[idSala] = (espOcupadoInf + salas->at(idSala).largura/2)*-1;
            espOcupadoInf += salas->at(idSala).largura;
        }
        if(espOcupadoSup < espOcupadoInf)
            espOcupadoMenor = espOcupadoSup;
        else
            espOcupadoMenor = espOcupadoInf;

        salasInseridas[idSala] = true;
        numSalasInseridas++;
    }

    delete idCandidatos;
    delete fluxoCandidatos;
    delete[] salasInseridas;
    return corredor;
}

double* gulosoRandomizado(float alpha, vector<Sala> *salas) {
    double menorCusto = INFINITY;
    double custo;
    double *menorCorredor = nullptr;
    double *corredor;
    for(int i=0; i<30; i++) {
        corredor = auxGulosoRandomizado(i, alpha, salas);
        custo = calculaCusto(salas, corredor);
        if(custo < menorCusto) {
            if(menorCorredor != nullptr)
                delete[] menorCorredor;
            menorCorredor = corredor;
            menorCusto = custo;
        } else
            delete[] corredor;
    }
    return menorCorredor;
}

double* gulosoReativo(vector<Sala> *salas, int numInteracoes, int tamanhoBloco) {
    vector<double> alphas;
    vector<double> probAlphas;
    vector<double> somatorioCustos;
    vector<int> quantSelecoes;
    double menorCusto = INFINITY;
    double custo;
    double *menorCorredor = nullptr;
    double *corredor;
    double auxTotalMedia;
    double acumulada;
    float melhorAlpha;
    int countBloco;
    int aleatorio;
    int indiceAlpha;

    for(int i=1; i<=NUMALPHAS; i++) {
        alphas.push_back(i*0.05);
        corredor = gulosoRandomizado(alphas[i-1], salas);
        custo = calculaCusto(salas, corredor);
        somatorioCustos.push_back(custo);
        quantSelecoes.push_back(1);
    }

    for(int i=0; i<NUMALPHAS; i++)
        probAlphas.push_back(100.0/NUMALPHAS);

    countBloco = 0;
    for(int i=NUMALPHAS; i<numInteracoes; i++) {
        if(countBloco == tamanhoBloco) {
            countBloco = 0;
            auxTotalMedia = 0;
            for(int j=0; j<alphas.size(); j++)
                auxTotalMedia += 1.0/(somatorioCustos[j]/quantSelecoes[j]);
            for(int j=0; j<alphas.size(); j++)
                probAlphas[j] = ((1.0/(somatorioCustos[j]/quantSelecoes[j])) / auxTotalMedia) * 100;
        }

        aleatorio = rand() % 100;
        acumulada = 0;
        for(int j=0; j<alphas.size(); j++) {
            acumulada += probAlphas[j];
            if(aleatorio < acumulada) {
                indiceAlpha = j;
                break;
            }
        }

        corredor = gulosoRandomizado(alphas[indiceAlpha], salas);
        custo = calculaCusto(salas, corredor);
        if(custo < menorCusto) {
            if(menorCorredor != nullptr)
                delete[] menorCorredor;
            menorCorredor = corredor;
            menorCusto = custo;
            melhorAlpha = alphas[indiceAlpha];
        } else
            delete[] corredor;

        indiceAlpha++;
    }
    cout << "Melhor Alpha: " << melhorAlpha << endl;
    return menorCorredor;
}

void cenarioUm(string arquivo) {
    double *solucao;
    vector<Sala> *salas;
    clock_t tempo[2];
    clock_t duracao;

    cout << "INSTANCIA: " << arquivo << endl << endl;

    salas = carregaInstancia("../Instancias/" + arquivo);
    tempo[0] = clock();
    solucao = guloso(salas);
    tempo[1] = clock();
    duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    cout << "GULOSO - Tempo: " << (float)duracao << "ms" << endl;
    apresentaSolucao(salas, solucao);
    cout << endl;
    delete[] solucao;

    tempo[0] = clock();
    solucao = gulosoRandomizado(0.3, salas);
    tempo[1] = clock();
    duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    cout << "GULOSO RANDOMIZADO - Tempo: " << (float)duracao << "ms" << endl;
    apresentaSolucao(salas, solucao);
    cout << endl;
    delete[] solucao;

    tempo[0] = clock();
    solucao = gulosoReativo(salas, 200, 20);
    tempo[1] = clock();
    duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    cout << "GULOSO REATIVO - Tempo: " << (float)duracao << "ms" << endl;
    apresentaSolucao(salas, solucao);
    cout << endl;
    delete[] solucao;

    for(int i=0; i<salas->size(); i++)
        delete[] salas->at(i).fluxo;
        
    delete salas;
}

int main()
{
    cenarioUm("Inst-10salas - 1374.txt");
    cenarioUm("Inst-11salas - 3439.txt");
    cenarioUm("Inst-56salas - 296220.txt");
    return 0;
}
