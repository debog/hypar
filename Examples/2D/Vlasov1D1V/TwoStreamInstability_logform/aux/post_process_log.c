/* This code extracts the maximum E-field for the Landau
 * damping case from the electric field files efield_<nnnnn>.dat,
 * which are text files with three columns (grid index, x-coordinate,
 * and the e-field value). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void IncrementFilename(char *f)
{
  if (f[7] == '9') {
    f[7] = '0';
    if (f[6] == '9') {
      f[6] = '0';
      if (f[5] == '9') {
        f[5] = '0';
        if (f[4] == '9') {
          f[4] = '0';
          if (f[3] == '9') {
            f[3] = '0';
            fprintf(stderr,"Warning: file increment hit max limit. Resetting to zero.\n");
          } else {
            f[3]++;
          }
        } else {
          f[4]++;
        }
      } else {
        f[5]++;
      }
    } else {
      f[6]++;
    }
  } else {
    f[7]++;
  }
}


int main() {
    int NI, NJ, ndims;
    double dt, samay;
    int file_op_iter, restart_iter = 0;
    char filename[50];

    FILE *in = fopen("solver.inp", "r");
    if (!in) {
        fprintf(stderr, "Error: File \"solver.inp\" not found.\n");
        return 1;
    } else {
        char word[100];
        fscanf(in, "%s", word);
        if (!strcmp(word, "begin")) {
            while (strcmp(word, "end")) {
                fscanf(in, "%s", word);
                if (!strcmp(word, "ndims")) fscanf(in, "%d", &ndims);
                else if (!strcmp(word, "size")) {
                    fscanf(in, "%d", &NI);
                    fscanf(in, "%d", &NJ);
                } else if (!strcmp(word, "dt")) fscanf(in, "%lf", &dt);
                else if (!strcmp(word, "file_op_iter")) fscanf(in, "%d", &file_op_iter);
                else if (!strcmp(word, "restart_iter")) fscanf(in, "%d", &restart_iter);
            }
        }
        fclose(in);
    }

    samay = (double)restart_iter * dt;
    strcpy(filename, "op_00000.dat");
    for (int t = 0; t < restart_iter; t++) {
        if ((t + 1) % file_op_iter == 0) IncrementFilename(filename);
    }

    while (1) {
        char out_filename[50];
        strcpy(out_filename, filename);
        out_filename[0] = 'p';  // Change the prefix from 'o' to 'p'

        FILE *in_file = fopen(filename, "r");
        if (!in_file) {
            fprintf(stderr, "Error: File \"%s\" not found.\n", filename);
            break;
        }

        FILE *out_file = fopen(out_filename, "w");
        if (!out_file) {
            fprintf(stderr, "Error: Could not open \"%s\" for writing.\n", out_filename);
            fclose(in_file);
            break;
        }

        char buffer[256];
        // Skip two header lines
        fgets(buffer, sizeof(buffer), in_file);
        fputs(buffer, out_file); // Copy the first header line
        fgets(buffer, sizeof(buffer), in_file);
        fputs(buffer, out_file); // Copy the second header line

        int i, j;
        double x, y, value;
        while (fscanf(in_file, "%d %d %lf %lf %lf", &i, &j, &x, &y, &value) == 5) {
            double exp_value = exp(value);
            fprintf(out_file, "%d %d %1.16e %1.16e %1.16e\n", i, j, x, y, exp_value);
        }

        fclose(in_file);
        fclose(out_file);

        IncrementFilename(filename);  // Increment to the next input and output file
    }

    return 0;
}
