      Save state for plotting in real time
        if (plotRealTime || t >= tEnd) {
            char filename[50];
            sprintf(filename, "%s/output_%d.txt",dirName, outputCount);
            save_to_file(U, mask, N, filename);
            outputCount++;
        }