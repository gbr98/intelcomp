#include <iostream>
#include <algorithm> 
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

struct AuxOrdena{
    int idCandidato;
    double fluxoCandidato;

    AuxOrdena(int idCandidato, double fluxoCandidato) {
        this->idCandidato = idCandidato;
        this->fluxoCandidato = fluxoCandidato;
    }
};

struct AuxOrdenaCorredor{
    double x;
    Sala* sala;

    AuxOrdenaCorredor(double x, Sala* sala) {
        this->x = x;
        this->sala = sala;
    }
};

vector<Sala>* carregaInstancia(string nomeArquivo) {
    vector<Sala> *salas;
    int numSalas;
    int *fluxo;
    ifstream arquivoEntrada;
    string str;
    stringstream ss;
    arquivoEntrada.open(nomeArquivo);

    if(arquivoEntrada.is_open()) {
        salas = new vector<Sala>;

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
    } else {
        cout << "Arquivo nao Encontrado" << endl;
        return nullptr;
    }
    return salas;
}

bool compara_sort(AuxOrdena* a, AuxOrdena* b) {
    return (a->fluxoCandidato < b->fluxoCandidato);
}

bool compara_sort_corredor(AuxOrdenaCorredor* a, AuxOrdenaCorredor* b) {
    return (a->x < b->x);
}

bool verificaSolucao(vector<Sala> *salas, double *corredor) {
    int numSalas = salas->size();
    int indiceSup;
    bool valida;
    double posicaoAtual;
    vector<AuxOrdenaCorredor*> auxCorredor;
    
    for(int i=0; i<numSalas; i++) 
        auxCorredor.push_back(new AuxOrdenaCorredor(corredor[i], &salas->at(i)));
    sort(auxCorredor.begin(), auxCorredor.end(), compara_sort_corredor);

    for(int i=0; i<numSalas; i++)
        if(auxCorredor[i]->x > 0) {
            indiceSup = i;
            break;
        }

    valida = true;

    //CORREDOR INFERIOR
    posicaoAtual = 0;
    for(int i=indiceSup-1; i>=0; i--) {
        if(posicaoAtual - (auxCorredor[i]->sala->largura / 2.0) != auxCorredor[i]->x) {
            valida = false;
            break;
        }
        posicaoAtual -= auxCorredor[i]->sala->largura;
    }

    //CORREDOR SUPERIOR
    if(valida) {
        posicaoAtual = 0;
        for(int i=indiceSup; i<numSalas; i++) {
            if(posicaoAtual + (auxCorredor[i]->sala->largura / 2.0) != auxCorredor[i]->x) {
                valida = false;
                break;
            }
            posicaoAtual += auxCorredor[i]->sala->largura;
        }
    }

    return valida;
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
    cout << "Custo: " << calculaCusto(salas, corredor);
    if(verificaSolucao(salas, corredor))
        cout << ", Valida" << endl;
    else
        cout << ", Invalida" << endl;
    cout << endl;
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

void buscaLocalSwap(vector<Sala> *salas, double* corredor) {
    int numSalas = salas->size();
    int indiceSup;
    bool atualizou;
    double diferenca;
    double custo;
    double novoCusto;
    vector<AuxOrdenaCorredor*> novoCorredor;
    vector<vector<double>*> fluxos;
    Sala* auxSwap;

    custo = calculaCusto(salas, corredor);

    for(int i=0; i<numSalas; i++) {
        novoCorredor.push_back(new AuxOrdenaCorredor(corredor[i], &salas->at(i)));
        fluxos.push_back(new vector<double>());
        for(int j=0; j<numSalas; j++) {
            fluxos[i]->push_back(salas->at(i).fluxo[j]);
        }
    }

    atualizou = true;
    while(atualizou) {
        atualizou = false;
        sort(novoCorredor.begin(), novoCorredor.end(), compara_sort_corredor);

        for(int i=0; i<numSalas; i++)
            if(novoCorredor[i]->x > 0) {
                indiceSup = i;
                break;
            }
        
        //SWAP INFERIOR COM INFERIOR
        for(int i=0; i<indiceSup; i++) {
            for(int j=i+1; j<indiceSup; j++) {
                diferenca = novoCorredor[i]->sala->largura - novoCorredor[j]->sala->largura;
                auxSwap = novoCorredor[i]->sala;
                novoCorredor[i]->sala = novoCorredor[j]->sala;
                novoCorredor[j]->sala = auxSwap;
            
                novoCorredor[j]->x -= diferenca/2;
                for(int k=j-1; k>i; k--) 
                    novoCorredor[k]->x -= diferenca;
                novoCorredor[i]->x -= diferenca/2;

                novoCusto = 0;
                for(int k=0; k<novoCorredor.size(); k++) {
                    for(int l=k+1; l<novoCorredor.size(); l++)
                            novoCusto += fabs(fabs(novoCorredor[k]->x) - fabs(novoCorredor[l]->x))*fluxos[novoCorredor[k]->sala->id]->at(novoCorredor[l]->sala->id);
                }
                
                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novoCorredor[i]->sala;
                    novoCorredor[i]->sala = novoCorredor[j]->sala;
                    novoCorredor[j]->sala = auxSwap;

                    diferenca *= -1;
                    novoCorredor[j]->x -= diferenca/2;
                    for(int k=j-1; k>i; k--) 
                        novoCorredor[k]->x -= diferenca;
                    novoCorredor[i]->x -= diferenca/2;

                    novoCusto = 0;
                    for(int k=0; k<novoCorredor.size(); k++) {
                        for(int l=k+1; l<novoCorredor.size(); l++)
                                novoCusto += fabs(fabs(novoCorredor[k]->x) - fabs(novoCorredor[l]->x))*fluxos[novoCorredor[k]->sala->id]->at(novoCorredor[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
        
        //SWAP SUPERIOR COM SUPERIOR
        for(int i=indiceSup; i<numSalas; i++) { 
            for(int j=i+1; j<numSalas; j++) {
                diferenca = novoCorredor[j]->sala->largura - novoCorredor[i]->sala->largura;
                auxSwap = novoCorredor[i]->sala;
                novoCorredor[i]->sala = novoCorredor[j]->sala;
                novoCorredor[j]->sala = auxSwap;
            
                novoCorredor[i]->x += diferenca/2;
                for(int k=i+1; k<j; k++)
                    novoCorredor[k]->x += diferenca;
                novoCorredor[j]->x += diferenca/2;

                novoCusto = 0;
                for(int k=0; k<novoCorredor.size(); k++) {
                    for(int l=k+1; l<novoCorredor.size(); l++)
                        novoCusto += fabs(fabs(novoCorredor[k]->x) - fabs(novoCorredor[l]->x))*fluxos[novoCorredor[k]->sala->id]->at(novoCorredor[l]->sala->id);
                }

                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novoCorredor[i]->sala;
                    novoCorredor[i]->sala = novoCorredor[j]->sala;
                    novoCorredor[j]->sala = auxSwap;
                
                    diferenca *= -1;
                    novoCorredor[i]->x += diferenca/2;
                    for(int k=i+1; k<j; k++)
                        novoCorredor[k]->x += diferenca;
                    novoCorredor[j]->x += diferenca/2;

                    novoCusto = 0;
                    for(int k=0; k<novoCorredor.size(); k++) {
                        for(int l=k+1; l<novoCorredor.size(); l++)
                            novoCusto += fabs(fabs(novoCorredor[k]->x) - fabs(novoCorredor[l]->x))*fluxos[novoCorredor[k]->sala->id]->at(novoCorredor[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
        
        //SWAP INFERIOR COM SUPERIOR
        for(int i=0; i<indiceSup; i++) {
            for(int j=indiceSup; j<numSalas; j++) {
                diferenca = novoCorredor[i]->sala->largura - novoCorredor[j]->sala->largura;
                auxSwap = novoCorredor[i]->sala;
                novoCorredor[i]->sala = novoCorredor[j]->sala;
                novoCorredor[j]->sala = auxSwap;
            
                novoCorredor[i]->x += diferenca/2;
                novoCorredor[j]->x += diferenca/2;
                for(int k=i-1; k>=0; k--) 
                    novoCorredor[k]->x += diferenca;
                for(int k=j+1; k<numSalas; k++) 
                    novoCorredor[k]->x += diferenca;


                novoCusto = 0;
                for(int k=0; k<novoCorredor.size(); k++) {
                    for(int l=k+1; l<novoCorredor.size(); l++)
                            novoCusto += fabs(fabs(novoCorredor[k]->x) - fabs(novoCorredor[l]->x))*fluxos[novoCorredor[k]->sala->id]->at(novoCorredor[l]->sala->id);
                }
                
                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novoCorredor[i]->sala;
                    novoCorredor[i]->sala = novoCorredor[j]->sala;
                    novoCorredor[j]->sala = auxSwap;
                
                    diferenca *= -1;
                    novoCorredor[i]->x += diferenca/2;
                    novoCorredor[j]->x += diferenca/2;
                    for(int k=i-1; k>=0; k--) 
                        novoCorredor[k]->x += diferenca;
                    for(int k=j+1; k<numSalas; k++) 
                        novoCorredor[k]->x += diferenca;

                    novoCusto = 0;
                    for(int k=0; k<novoCorredor.size(); k++) {
                        for(int l=k+1; l<novoCorredor.size(); l++)
                                novoCusto += fabs(fabs(novoCorredor[k]->x) - fabs(novoCorredor[l]->x))*fluxos[novoCorredor[k]->sala->id]->at(novoCorredor[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
    }

    for(int i=0; i<numSalas; i++)
        corredor[novoCorredor[i]->sala->id] = novoCorredor[i]->x;

    for(int i=0; i<numSalas; i++) {
        delete fluxos[i];
        delete novoCorredor[i];
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
        maiorFluxo = -1;
        for(int i=0; i<numSalas; i++) {
            if(!salasInseridas[i]) {
                custoAux = fabs((fabs(corredor[idUltimaAdd]) - espOcupadoMenor + salas->at(i).largura/2.0)) * salas->at(idUltimaAdd).fluxo[i];
                if(custoAux > maiorFluxo) {
                    idMaiorFluxo = i;
                    maiorFluxo = custoAux;
                    larguraMaiorFluxo = salas->at(i).largura;
                }
            }
        }
        
        if(espOcupadoSup < espOcupadoInf) {
            corredor[idMaiorFluxo] = espOcupadoSup + larguraMaiorFluxo/2.0;
            espOcupadoSup += larguraMaiorFluxo;
        } else {
            corredor[idMaiorFluxo] = (espOcupadoInf + larguraMaiorFluxo/2.0)*-1;
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
    int aleatorio;
    int idSala;
    int idUltimaAdd;
    double *corredor = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    vector<AuxOrdena*> candidatos;

    for(int i=0; i<numSalas; i++) {
        corredor[i] = 0;
        candidatos.push_back(new AuxOrdena(i, 0));
    }

    idSala = selecionaOrigem(salas);
    espOcupadoSup = salas->at(idSala).largura;
    corredor[idSala] = espOcupadoSup/2.0;
    candidatos.erase(candidatos.begin()+idSala);

    idUltimaAdd = idSala;
    while(candidatos.size() > 0) {
        for(int i=0; i<numSalas; i++)
            candidatos[i]->fluxoCandidato = fabs((fabs(corredor[idUltimaAdd]) - espOcupadoMenor + salas->at(candidatos[i]->idCandidato).largura/2.0)) * salas->at(idUltimaAdd).fluxo[candidatos[i]->idCandidato];
        
        sort(candidatos.begin(), candidatos.end(), compara_sort);
        
        aleatorio = rand()%(int)ceil(alpha*candidatos.size());
        idSala = candidatos[aleatorio]->idCandidato;
        if(espOcupadoSup < espOcupadoInf) {
            corredor[idSala] = espOcupadoSup + salas->at(idSala).largura/2.0;
            espOcupadoSup += salas->at(idSala).largura;
        } else {
            corredor[idSala] = (espOcupadoInf + salas->at(idSala).largura/2.0)*-1.0;
            espOcupadoInf += salas->at(idSala).largura;
        }
        if(espOcupadoSup < espOcupadoInf)
            espOcupadoMenor = espOcupadoSup;
        else
            espOcupadoMenor = espOcupadoInf;

        delete candidatos[aleatorio];
        candidatos.erase(candidatos.begin()+aleatorio);
    }

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
    return menorCorredor;
}

double* solucaoInicial(vector<Sala> *salas) {
    double *solucao;
    double *melhorSolucao;
    double menorCusto;
    double custo;
    clock_t tempo[2];
    clock_t duracao;

    tempo[0] = clock();
    solucao = guloso(salas);
    tempo[1] = clock();
    duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    //cout << "GULOSO - Tempo: " << (float)duracao << "ms" << endl;
    //apresentaSolucao(salas, solucao);
    buscaLocalSwap(salas, solucao);
    //apresentaSolucao(salas, solucao);

    melhorSolucao = solucao;
    menorCusto = calculaCusto(salas, solucao);

    tempo[0] = clock();
    solucao = gulosoRandomizado(0.3, salas);
    tempo[1] = clock();
    duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    //cout << "GULOSO RANDOMIZADO - Tempo: " << (float)duracao << "ms" << endl;
    //apresentaSolucao(salas, solucao);
    buscaLocalSwap(salas, solucao);
    //apresentaSolucao(salas, solucao);
    
    custo = calculaCusto(salas, solucao);
    if(custo < menorCusto) {
        menorCusto = custo;
        delete melhorSolucao;
        melhorSolucao = solucao;
    } else
        delete solucao;

    tempo[0] = clock();
    solucao = gulosoReativo(salas, 200, 20);
    tempo[1] = clock();
    duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    //cout << "GULOSO REATIVO - Tempo: " << (float)duracao << "ms" << endl;
    //apresentaSolucao(salas, solucao);
    buscaLocalSwap(salas, solucao);
    //apresentaSolucao(salas, solucao);
    
    custo = calculaCusto(salas, solucao);
    if(custo < menorCusto) {
        menorCusto = custo;
        delete melhorSolucao;
        melhorSolucao = solucao;
    } else
        delete solucao;

    return melhorSolucao;
}

void cenarioUm(string arquivo, double menorCusto) {
    vector<Sala> *salas;
    double* solucao;
    double custo;

    cout << "INSTANCIA: " << arquivo << endl;
    salas = carregaInstancia("../Instancias/" + arquivo);
    if(salas == nullptr)
        return;

    solucao = solucaoInicial(salas);
    custo = calculaCusto(salas, solucao);

    cout << "Melhor Solucao: ";
    for(int i=0; i<salas->size(); i++)
        cout << solucao[i] << " ";
    cout << endl << "Custo: " << custo << ", Erro: " << (custo - menorCusto)*100/menorCusto << endl << endl;

    for(int i=0; i<salas->size(); i++)
        delete[] salas->at(i).fluxo;

    delete solucao;
    delete salas;
}

int main()
{
    cenarioUm("Inst-10salas - 1374.txt", 1374);
    cenarioUm("Inst-11salas - 3439.txt", 3439);
    cenarioUm("Inst-56salas - 296220.txt", 296220);

    return 0;
}
