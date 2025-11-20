#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <string.h>
#include <cstdlib>
#include <dirent.h>

#include <iostream>

#include "../include/define.h"
#include "../include/evolve.h"
#include "../include/io.h"
#include "../include/monte_carlo.h"
#include "../include/parameter_order.h"
#include "../include/parameters.h"
#include "../include/potential.h"
#include "../include/simulator.h"
#define MCLICS_VERSION "0.1"

void cleanPreviousOutput() {
    printf("Limpando arquivos de simulações anteriores...\n");
    
    // Remove arquivos específicos
    const char* specificFiles[] = {"ic.csv", "po.dat", NULL};
    for (int i = 0; specificFiles[i] != NULL; i++) {
        remove(specificFiles[i]);
        printf("Removido: %s\n", specificFiles[i]);
    }
    
    // Remove arquivos director_field_*.csv manualmente
    DIR *dir = opendir(".");
    if (dir != NULL) {
        struct dirent *entry;
        while ((entry = readdir(dir)) != NULL) {
            std::string filename = entry->d_name;
            // Verifica se o arquivo começa com "director_field_" e termina com ".csv"
            if (filename.find("director_field_") == 0 && 
                filename.length() >= 4 &&
                filename.substr(filename.length() - 4) == ".csv") {
                remove(entry->d_name);
                printf("Removido: %s\n", entry->d_name);
            }
        }
        closedir(dir);
    }
    
    printf("Limpeza concluída!\n\n");
}

void printSimulationSummary(Parameters *params) {
    printf("\n=== RELATÓRIO FINAL DA SIMULAÇÃO ===\n");
    
    FILE *poFile = fopen("po.dat", "r");
    if (poFile) {
        double sumS = 0.0, sumE = 0.0;
        double maxS = -1e9, minS = 1e9;
        double maxE = -1e9, minE = 1e9;
        int count = 0;
        char line[256];
        
        // Pular cabeçalho (verificando retorno)
        if (fgets(line, sizeof(line), poFile) == NULL) {
            // Arquivo vazio ou erro
            fclose(poFile);
            printf("Arquivo po.dat vazio ou erro na leitura\n");
            return;
        }
        
        // Ler dados
        double T, S, varS, E, varE;
        while (fgets(line, sizeof(line), poFile) != NULL) {
            if (sscanf(line, "%lf %lf %lf %lf %lf", &T, &S, &varS, &E, &varE) == 5) {
                sumS += S;
                sumE += E;
                if (S > maxS) maxS = S;
                if (S < minS) minS = S;
                if (E > maxE) maxE = E;
                if (E < minE) minE = E;
                count++;
            }
        }
        fclose(poFile);
        
        if (count > 0) {
            printf("Parâmetro de Ordem (S):\n");
            printf("  Média: %.4f\n", sumS / count);
            printf("  Máximo: %.4f\n", maxS);
            printf("  Mínimo: %.4f\n", minS);
            
            printf("\nEnergia (E):\n");
            printf("  Média: %.4f\n", sumE / count);
            printf("  Máxima: %.4f\n", maxE);
            printf("  Mínima: %.4f\n", minE);
            
            printf("\nEstatísticas:\n");
            printf("  Número de temperaturas: %d\n", count);
        } else {
            printf("Nenhum dado válido encontrado em po.dat\n");
        }
    } else {
        printf("Arquivo po.dat não encontrado\n");
    }
    
    printf("\nConfiguração do Sistema:\n");
    printf("  Dimensões: %dx%dx%d\n", params->Nx, params->Ny, params->Nz);
    printf("  Potencial: %s\n", params->potential);
    printf("  Geometria: %s\n", params->geometry);
    printf("  Condições de contorno: %s/%s/%s\n", 
           params->XBoundtype, params->YBoundtype, params->ZBoundtype);
    
    printf("====================================\n\n");
}

int main(int argc, char **argv) {
  printf("### Starting McLiCS version: %s ###\n\n",MCLICS_VERSION);
  cleanPreviousOutput();
  Parameters params = read_input_file(argv[1]);
  print_parameters(params);
  char fname[1000];
  simulator *sim = new simulator(&params);
  sim->Setup_simmulation(params);

  sprintf(fname, "ic.csv");
  sim->print_n(fname, &params);
  sim->evolve->run();

  printSimulationSummary(&params);
  delete sim;

  return 0;
}
