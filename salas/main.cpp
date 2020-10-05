#include <iostream>
#include <algorithm> 
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <time.h>

#define NUMALPHAS 3

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

struct AuxOrdenaSolucao{
    double x;
    Sala* sala;

    AuxOrdenaSolucao(double x, Sala* sala) {
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

bool compara_sort_solucao(AuxOrdenaSolucao* a, AuxOrdenaSolucao* b) {
    return (a->x < b->x);
}

bool verificaSolucao(vector<Sala> *salas, double *solucao) {
    int numSalas = salas->size();
    int indiceSup;
    bool valida;
    double posicaoAtual;
    vector<AuxOrdenaSolucao*> auxSolucao;
    
    for(int i=0; i<numSalas; i++) 
        auxSolucao.push_back(new AuxOrdenaSolucao(solucao[i], &salas->at(i)));
    sort(auxSolucao.begin(), auxSolucao.end(), compara_sort_solucao);

    for(int i=0; i<numSalas; i++)
        if(auxSolucao[i]->x > 0) {
            indiceSup = i;
            break;
        }

    valida = true;

    //CORREDOR INFERIOR
    posicaoAtual = 0;
    for(int i=indiceSup-1; i>=0; i--) {
        if(posicaoAtual - (auxSolucao[i]->sala->largura / 2.0) != auxSolucao[i]->x) {
            valida = false;
            break;
        }
        posicaoAtual -= auxSolucao[i]->sala->largura;
    }

    //CORREDOR SUPERIOR
    if(valida) {
        posicaoAtual = 0;
        for(int i=indiceSup; i<numSalas; i++) {
            if(posicaoAtual + (auxSolucao[i]->sala->largura / 2.0) != auxSolucao[i]->x) {
                valida = false;
                break;
            }
            posicaoAtual += auxSolucao[i]->sala->largura;
        }
    }

    return valida;
}

double calculaCusto(vector<Sala> *salas, double *solucao) {
    double custo = 0;

    if(solucao == nullptr)
        return -1;

    for(int i=0; i<salas->size(); i++) {
        for(int j=i+1; j<salas->size(); j++)
                custo += fabs(fabs(solucao[i]) - fabs(solucao[j]))*salas->at(i).fluxo[j];
    }
    return custo;
}

double estimaCusto(vector<Sala> *salas, double custo, double *solucao, double posCandidato, int idCandidato) {
    for(int i=0; i<salas->at(0).numSalas; i++) {
        if(solucao[i] != 0)
            custo += fabs(fabs(solucao[i]) - fabs(posCandidato))*salas->at(i).fluxo[idCandidato];
    }
    return custo;
}

void apresentaSolucao(vector<Sala> *salas, double* solucao) {
    cout << "Custo: " << calculaCusto(salas, solucao);
    if(verificaSolucao(salas, solucao))
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

void buscaLocalSwap(vector<Sala> *salas, double* solucao) {
    int numSalas = salas->size();
    int indiceSup;
    bool atualizou;
    double diferenca;
    double custo;
    double novoCusto;
    vector<AuxOrdenaSolucao*> novaSolucao;
    vector<vector<double>*> fluxos;
    Sala* auxSwap;

    custo = calculaCusto(salas, solucao);

    for(int i=0; i<numSalas; i++) {
        novaSolucao.push_back(new AuxOrdenaSolucao(solucao[i], &salas->at(i)));
        fluxos.push_back(new vector<double>());
        for(int j=0; j<numSalas; j++) {
            fluxos[i]->push_back(salas->at(i).fluxo[j]);
        }
    }

    atualizou = true;
    while(atualizou) {
        atualizou = false;
        sort(novaSolucao.begin(), novaSolucao.end(), compara_sort_solucao);

        for(int i=0; i<numSalas; i++)
            if(novaSolucao[i]->x > 0) {
                indiceSup = i;
                break;
            }
        
        //SWAP INFERIOR COM INFERIOR
        for(int i=0; i<indiceSup; i++) {
            for(int j=i+1; j<indiceSup; j++) {
                diferenca = novaSolucao[i]->sala->largura - novaSolucao[j]->sala->largura;
                auxSwap = novaSolucao[i]->sala;
                novaSolucao[i]->sala = novaSolucao[j]->sala;
                novaSolucao[j]->sala = auxSwap;
            
                novaSolucao[j]->x -= diferenca/2;
                for(int k=j-1; k>i; k--) 
                    novaSolucao[k]->x -= diferenca;
                novaSolucao[i]->x -= diferenca/2;

                novoCusto = 0;
                for(int k=0; k<novaSolucao.size(); k++) {
                    for(int l=k+1; l<novaSolucao.size(); l++)
                            novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                }
                
                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novaSolucao[i]->sala;
                    novaSolucao[i]->sala = novaSolucao[j]->sala;
                    novaSolucao[j]->sala = auxSwap;

                    diferenca *= -1;
                    novaSolucao[j]->x -= diferenca/2;
                    for(int k=j-1; k>i; k--) 
                        novaSolucao[k]->x -= diferenca;
                    novaSolucao[i]->x -= diferenca/2;

                    novoCusto = 0;
                    for(int k=0; k<novaSolucao.size(); k++) {
                        for(int l=k+1; l<novaSolucao.size(); l++)
                                novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
        
        //SWAP SUPERIOR COM SUPERIOR
        for(int i=indiceSup; i<numSalas; i++) { 
            for(int j=i+1; j<numSalas; j++) {
                diferenca = novaSolucao[j]->sala->largura - novaSolucao[i]->sala->largura;
                auxSwap = novaSolucao[i]->sala;
                novaSolucao[i]->sala = novaSolucao[j]->sala;
                novaSolucao[j]->sala = auxSwap;
            
                novaSolucao[i]->x += diferenca/2;
                for(int k=i+1; k<j; k++)
                    novaSolucao[k]->x += diferenca;
                novaSolucao[j]->x += diferenca/2;

                novoCusto = 0;
                for(int k=0; k<novaSolucao.size(); k++) {
                    for(int l=k+1; l<novaSolucao.size(); l++)
                        novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                }

                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novaSolucao[i]->sala;
                    novaSolucao[i]->sala = novaSolucao[j]->sala;
                    novaSolucao[j]->sala = auxSwap;
                
                    diferenca *= -1;
                    novaSolucao[i]->x += diferenca/2;
                    for(int k=i+1; k<j; k++)
                        novaSolucao[k]->x += diferenca;
                    novaSolucao[j]->x += diferenca/2;

                    novoCusto = 0;
                    for(int k=0; k<novaSolucao.size(); k++) {
                        for(int l=k+1; l<novaSolucao.size(); l++)
                            novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
        
        //SWAP INFERIOR COM SUPERIOR
        for(int i=0; i<indiceSup; i++) {
            for(int j=indiceSup; j<numSalas; j++) {
                diferenca = novaSolucao[i]->sala->largura - novaSolucao[j]->sala->largura;
                auxSwap = novaSolucao[i]->sala;
                novaSolucao[i]->sala = novaSolucao[j]->sala;
                novaSolucao[j]->sala = auxSwap;
            
                novaSolucao[i]->x += diferenca/2;
                novaSolucao[j]->x += diferenca/2;
                for(int k=i-1; k>=0; k--) 
                    novaSolucao[k]->x += diferenca;
                for(int k=j+1; k<numSalas; k++) 
                    novaSolucao[k]->x += diferenca;


                novoCusto = 0;
                for(int k=0; k<novaSolucao.size(); k++) {
                    for(int l=k+1; l<novaSolucao.size(); l++)
                            novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                }
                
                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novaSolucao[i]->sala;
                    novaSolucao[i]->sala = novaSolucao[j]->sala;
                    novaSolucao[j]->sala = auxSwap;
                
                    diferenca *= -1;
                    novaSolucao[i]->x += diferenca/2;
                    novaSolucao[j]->x += diferenca/2;
                    for(int k=i-1; k>=0; k--) 
                        novaSolucao[k]->x += diferenca;
                    for(int k=j+1; k<numSalas; k++) 
                        novaSolucao[k]->x += diferenca;

                    novoCusto = 0;
                    for(int k=0; k<novaSolucao.size(); k++) {
                        for(int l=k+1; l<novaSolucao.size(); l++)
                                novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
    }

    for(int i=0; i<numSalas; i++)
        solucao[novaSolucao[i]->sala->id] = novaSolucao[i]->x;

    for(int i=0; i<numSalas; i++) {
        delete fluxos[i];
        delete novaSolucao[i];
    }
}

double* guloso(vector<Sala> *salas) {
    int numSalas = salas->size();
    int numSalasInseridas = 0;
    int idAtual;
    int idMaiorFluxo;
    int idUltimaAdd;
    double maiorFluxo;
    double *solucao = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    double larguraMaiorFluxo;
    double custo = 0;
    double custoAux;
    bool *salasInseridas = new bool[numSalas];

    for(int i=0; i<numSalas; i++) {
        solucao[i] = 0;
        salasInseridas[i] = false;
    }

    idAtual = selecionaOrigem(salas);
    espOcupadoSup = salas->at(idAtual).largura;
    solucao[idAtual] = espOcupadoSup/2;
    salasInseridas[idAtual] = true;
    numSalasInseridas++;

    idUltimaAdd = idAtual;
    while(numSalasInseridas < numSalas) {
        maiorFluxo = -1;
        for(int i=0; i<numSalas; i++) {
            if(!salasInseridas[i]) {
                custoAux = fabs(fabs(solucao[idUltimaAdd]) - (espOcupadoMenor + salas->at(i).largura/2.0)) * salas->at(idUltimaAdd).fluxo[i];
                if(custoAux > maiorFluxo) {
                    idMaiorFluxo = i;
                    maiorFluxo = custoAux;
                    larguraMaiorFluxo = salas->at(i).largura;
                }
            }
        }
        
        if(espOcupadoSup < espOcupadoInf) {
            solucao[idMaiorFluxo] = espOcupadoSup + larguraMaiorFluxo/2.0;
            espOcupadoSup += larguraMaiorFluxo;
        } else {
            solucao[idMaiorFluxo] = (espOcupadoInf + larguraMaiorFluxo/2.0)*-1;
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
    return solucao;
}

double* auxGulosoRandomizado(vector<Sala> *salas, float alpha, int seed) {
    srand(seed);
    int numSalas = salas->size();
    int aleatorio;
    int idSala;
    int idUltimaAdd;
    double *solucao = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    vector<AuxOrdena*> candidatos;

    for(int i=0; i<numSalas; i++) {
        solucao[i] = 0;
        candidatos.push_back(new AuxOrdena(i, 0));
    }

    idSala = selecionaOrigem(salas);
    espOcupadoSup = salas->at(idSala).largura;
    solucao[idSala] = espOcupadoSup/2.0;
    candidatos.erase(candidatos.begin()+idSala);

    idUltimaAdd = idSala;
    while(candidatos.size() > 0) {
        for(int i=0; i<candidatos.size(); i++)
            candidatos[i]->fluxoCandidato = fabs(fabs(solucao[idUltimaAdd]) - (espOcupadoMenor + salas->at(candidatos[i]->idCandidato).largura/2.0)) * salas->at(idUltimaAdd).fluxo[candidatos[i]->idCandidato];
        sort(candidatos.begin(), candidatos.end(), compara_sort);
        
        aleatorio = rand()%(int)ceil(alpha*candidatos.size());
        idSala = candidatos[aleatorio]->idCandidato;
        
        if(espOcupadoSup < espOcupadoInf) {
            solucao[idSala] = espOcupadoSup + salas->at(idSala).largura/2.0;
            espOcupadoSup += salas->at(idSala).largura;
        } else {
            solucao[idSala] = (espOcupadoInf + salas->at(idSala).largura/2.0)*-1.0;
            espOcupadoInf += salas->at(idSala).largura;
        }
        if(espOcupadoSup < espOcupadoInf)
            espOcupadoMenor = espOcupadoSup;
        else
            espOcupadoMenor = espOcupadoInf;

        delete candidatos[aleatorio];
        candidatos.erase(candidatos.begin()+aleatorio);
    }
    
    return solucao;
}

double* gulosoRandomizado(vector<Sala> *salas, float alpha, int seed) {
    double menorCusto = INFINITY;
    double custo;
    double *melhorSolucao = nullptr;
    double *solucao;
    for(int i=0; i<100; i++) {
        solucao = auxGulosoRandomizado(salas, alpha, seed+i);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
            menorCusto = custo;
        } else
            delete[] solucao;
    }
    return melhorSolucao;
}

double* gulosoReativo(vector<Sala> *salas, int numInteracoes, int tamanhoBloco) {
    vector<double> alphas;
    vector<double> probAlphas;
    vector<double> somatorioCustos;
    vector<int> quantSelecoes;
    double menorCusto = INFINITY;
    double custo;
    double *melhorSolucao = nullptr;
    double *solucao;
    double auxTotalMedia;
    double acumulada;
    float melhorAlpha;
    int countBloco;
    int aleatorio;
    int indiceAlpha;

    alphas.push_back(0.5);
    alphas.push_back(0.25);
    alphas.push_back(0.3);

    for(int i=0; i<NUMALPHAS; i++) {
        solucao = gulosoRandomizado(salas, alphas[i], 0);
        custo = calculaCusto(salas, solucao);
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

        solucao = gulosoRandomizado(salas, alphas[indiceAlpha], i);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
            menorCusto = custo;
            melhorAlpha = alphas[indiceAlpha];
        } else
            delete[] solucao;

        indiceAlpha++;
    }
    cout << "MELHOR: " << melhorAlpha << endl;
    return melhorSolucao;
}

double* graspGuloso(vector<Sala> *salas, int numInteracoes) {
    double *solucao;
    double *melhorSolucao;
    double menorCusto;
    double custo;

    melhorSolucao = nullptr;
    menorCusto = INFINITY;
    for(int i=0; i<numInteracoes; i++) {
        solucao = guloso(salas);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            menorCusto = custo;
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
        } else
            delete[] solucao;
    }

    return melhorSolucao;
}

double* graspRandomizado(vector<Sala> *salas, int numInteracoes, float alpha) {
    double *solucao;
    double *melhorSolucao;
    double menorCusto;
    double custo;

    melhorSolucao = nullptr;
    menorCusto = INFINITY;
    for(int i=0; i<numInteracoes; i++) {
        solucao = gulosoRandomizado(salas, alpha, i);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            menorCusto = custo;
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
        } else
            delete[] solucao;
    }

    return melhorSolucao;
}

double* graspReativo(vector<Sala> *salas, int numInteracoes, int numIteracoesReativo, int tamanhoBloco) {
    double *solucao;
    double *melhorSolucao;
    double menorCusto;
    double custo;

    melhorSolucao = nullptr;
    menorCusto = INFINITY;
    for(int i=0; i<numInteracoes; i++) {
        solucao = gulosoReativo(salas, numIteracoesReativo, tamanhoBloco);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            menorCusto = custo;
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
        } else
            delete[] solucao;
    }

    return melhorSolucao;
}

void cenarioUm(string arquivo, double custoObjetivo) {
    int melhorMetodo;
    double *solucao;
    double menorCusto;
    double custo[6];
    double tempoMetodos[6];
    clock_t tempo[2];
    string metodos[6] = {"Guloso", "Guloso Randomizado", "Guloso Reativo", "GRASP-Guloso", "GRASP-Randomizado", "GRASP-Reativo"};
    vector<Sala> *salas;

    salas = carregaInstancia("../Instancias/" + arquivo);
    if(salas == nullptr) 
        return;

    tempo[0] = clock();
    solucao = guloso(salas);
    tempo[1] = clock();
    tempoMetodos[0] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    custo[0] = calculaCusto(salas, solucao);
    menorCusto = custo[0];
    melhorMetodo = 0;
    delete[] solucao;
    
    tempo[0] = clock();
    solucao = gulosoRandomizado(salas, 0.1, 0);
    tempo[1] = clock();
    tempoMetodos[1] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    custo[1] = calculaCusto(salas, solucao);
    if(custo[1] < menorCusto) {
        menorCusto = custo[1];
        melhorMetodo = 1;
    }
    delete[] solucao;
    
    //O REATIVO DEMORA MUITO, ACHO MELHOR NÃO SUAR. EU USEI MAIS PRA AVALIAR OS ALPHAS MESMO
    /*
    tempo[0] = clock();
    solucao = gulosoReativo(salas, 200, 20);
    tempo[1] = clock();
    tempoMetodos[2] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    custo[2] = calculaCusto(salas, solucao);
    if(custo[2] < menorCusto) {
        menorCusto = custo[2];
        melhorMetodo = 2;
    }
    delete[] solucao;
    */

    tempoMetodos[2] = 0;
    custo[2] = 787878;

    tempo[0] = clock();
    solucao = graspGuloso(salas, 1);
    tempo[1] = clock();
    tempoMetodos[3] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    custo[3] = calculaCusto(salas, solucao);
    if(custo[3] < menorCusto) {
        menorCusto = custo[3];
        melhorMetodo = 3;
    }
    delete[] solucao;

    tempo[0] = clock();
    solucao = graspRandomizado(salas, 500, 0.7);
    //FUI MUDANDO O NÚMERO DE ITERAÇÕES DO GRASP E O ALPHA UTILIZADO. RESULTADOS ABAIXO PARA AS DUAS PRIMEIRAS INSTÂNCIAS
    //100 iterações
    //0.5 1393.5 3513.5
    //0.15 1434.5 3466.5 
    //0.25 1380.5 3548.5
    //0.3 1380.5 3545.5
    //0.35 1424.5 3466.5
    //0.4 1405.5 3513
    //0.55 1380.5 3504.5
    //0.6 1380.5 3535.5
    //0.65 1378.5 3580.5
    //0.85 1378.5 3606.5

    //200 iterações
    //0.5 1393.5 3477.5
    //0.3 1380.5 3439.5

    //500 iterações
    //0.65 1374.5 3466.5
    
    tempo[1] = clock();
    tempoMetodos[4] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    custo[4] = calculaCusto(salas, solucao);
    if(custo[4] < menorCusto) {
        menorCusto = custo[4];
        melhorMetodo = 4;
    }
    delete[] solucao;
    
    //DEMORA MUUUUUUITO
    /*
    tempo[0] = clock();
    solucao = graspReativo(salas, 10, 200, 50);
    tempo[1] = clock();
    tempoMetodos[5] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
    custo[5] = calculaCusto(salas, solucao);
    if(custo[5] < menorCusto) {
        menorCusto = custo[5];
        melhorMetodo = 5;
    }
    delete[] solucao;
    */
    tempoMetodos[5] = 0;
    custo[5] = 787878;

    for(int i=0; i<6; i++)
        cout << metodos[i] << ": " << "Custo: " << custo[i] << ", Erro: " << (custo[i] - custoObjetivo)*100/custoObjetivo << ", Tempo: " << tempoMetodos[i] << "ms" << endl;
    cout << "Melhor Metodo: " << metodos[melhorMetodo] << endl << endl;
    
    for(int i=0; i<salas->size(); i++)
        delete[] salas->at(i).fluxo;

    delete salas;
}

int main() {
    string arquivo;

    arquivo = "Inst-10salas - 1374.txt";
    cout << "INSTANCIA: " << arquivo << endl;
    cenarioUm(arquivo, 1374);

    arquivo = "Inst-11salas - 3439.txt";
    cout << "INSTANCIA: " << arquivo << endl;
    cenarioUm(arquivo, 3439);

    //DEIXEI ESSA INSTÂNCIA COMENTADA PQ ELA DEMORA MUITO PRA RODAR O GRASP. 56 SALAS PORRA
    /*
    arquivo = "Inst-56salas - 296220.txt";
    cout << "INSTANCIA: " << arquivo << endl;
    cenarioUm(arquivo, 296220);
    */
    return 0;
}