#include <stdio.h>

int main() {
    FILE *inputFile = fopen("solapamiento.dat", "r");
    FILE *outputFile = fopen("solapamientobien.dat", "w");

    int currentLabel = -1;
    double currentValue = 0.0;

    while (!feof(inputFile)) {
        int label;
        double value;

        if (fscanf(inputFile, " %d,%lf", &label, &value) != 2) {
            break;
        }

        if (currentLabel != label) {
            if (currentLabel != -1) {
                fprintf(outputFile, "\n");
            }
            currentLabel = label;
            currentValue = value;
        }

        fprintf(outputFile, " %d,%f\n", label, value);
    }

    rewind(inputFile);
    currentLabel = -1;

    while (!feof(inputFile)) {
        int label;
        double value;

        if (fscanf(inputFile, " %d,%lf", &label, &value) != 2) {
            break;
        }

        if (currentLabel != label) {
            if (currentLabel != -1) {
                fprintf(outputFile, "\n");
            }
            currentLabel = label;
            currentValue = value;
        }

        fprintf(outputFile, " %d,%f\n", label, -value);
    }

    fclose(inputFile);
    fclose(outputFile);

    return 0;
}