import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.*;
class DataAndSize {
    private double[][] data;
    private int numRows;
    private int numCols;

    public DataAndSize(double[][] data, int numRows, int numCols) {
        this.data = data;
        this.numRows = numRows;
        this.numCols = numCols;
    }

    public double[][] getData() {
        return data;
    }

    public int getNumRows() {
        return numRows;
    }

    public int getNumCols() {
        return numCols;
    }
}


public class grey_Wolf_Optimization {
    double r1;
    double r2;
    int N;
    int D;
    int maxiter;
    double[] alfa;
    double beta[];
    double delta[];
    double Lower[];
    double Upper[];

    double XX[][];
    double X1[][];
    double X2[][];
    double X3[][];
    double BESTVAL;
    double a[];
    double A1[];
    double C1[];
    double A2[];
    double C2[];
    double A3[];
    double C3[];

    public grey_Wolf_Optimization (int imaxiter) {// constractor
        maxiter=imaxiter;


        N=20;
        D=10;

        a=new double[D];
        XX=new double[N][D];
        alfa=new double[D];
        beta=new double[D];
        delta=new double[D];
        A1=new double[D];
        C1=new double[D];
        A2=new double[D];
        C2=new double[D];
        A3=new double[D];
        C3=new double[D];
        X1=new double[N][D];
        X2=new double[N][D];
        X3=new double[N][D];

    }

    double[][] sort_and_index(double[][] XXX) {


        String filePath = "/home/skstar/Desktop/HuntingTon_500*50.xlsx";
        DataAndSize dataAndSize = mainData(filePath);

        double[][] data = dataAndSize.getData();
        int numRows = dataAndSize.getNumRows();
        int numCols = dataAndSize.getNumCols();

        int[][] subsetIndices = selectSubsetIndices(numRows,numCols);
        double[] residueScores = calculateMeanSquaredResidue(data, subsetIndices);

        int N = subsetIndices.length;
        int D = subsetIndices[0].length;

        double[] fitnesValue=new double[residueScores.length];
        for(int i=0;i<residueScores.length;i++) {
            fitnesValue[i]=residueScores[i];
        }


        ArrayList<Double> newFitnessValue=new ArrayList<Double>(); //Copy the fitnessvaule on newFitnessValue
        for(int i=0;i<N;i++) {
            newFitnessValue.add(fitnesValue[i]);
        }

        ArrayList<Double> storeNewFitVal = new ArrayList<Double>(newFitnessValue);  // store newFitnessValue for latter use

        Collections.sort(newFitnessValue);

        double[] sortedValue=new double[newFitnessValue.size()];

        Iterator<Double> iterator=newFitnessValue.iterator();   //
        int ii=0;
        while(iterator.hasNext()) {
            sortedValue[ii]=iterator.next().doubleValue();
            ii++;
        }


        int[] indexes=new int[newFitnessValue.size()];
        for(int n=0;n<newFitnessValue.size(); n++) {
            indexes[n]=storeNewFitVal.indexOf(newFitnessValue.get(n));
        }
        double[][] B=new double[N][D];
        for(int i=0;i<N;i++) {
            for(int j=0;j<D;j++) {
                B[i][j]=XXX[indexes[i]][j];
            }
        }

        return B ;
    }







    DataAndSize mainData(String filePath) {
        double[][] data = null;
        int numRows = 0;
        int numCols = 0;

        try {
            // Load the Excel workbook
            File file = new File(filePath);
            InputStream inputStream = new FileInputStream(file);
            Workbook workbook = WorkbookFactory.create(inputStream);

            // Select the worksheet you want to read data from
            Sheet sheet = workbook.getSheetAt(0);

            // Retrieve the number of rows and columns
            numRows = sheet.getLastRowNum() + 1;
            numCols = sheet.getRow(0).getLastCellNum();

            // Initialize a 2D array with the same dimensions as the worksheet
            data = new double[numRows][numCols];

            // Iterate over each row and column in the worksheet and store the data in the 2D array
            for (int i = 0; i < numRows; i++) {
                Row row = sheet.getRow(i);
                for (int j = 0; j < numCols; j++) {
                    Cell cell = row.getCell(j);
                    double value = cell.getNumericCellValue();
                    data[i][j] = value;
                }
            }

            // Close the workbook and input stream
            workbook.close();
            inputStream.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

        return new DataAndSize(data, numRows, numCols);
    }













// ...

    int[][] selectSubsetIndices(int numRows, int numCols) {
        int[][] subsetIndices = new int[20][10];
        Random rand = new Random(1234); // seed value is 1234
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 10; j++) {
                if (j < 6) {
                    int randomRow = rand.nextInt(numRows);
                    subsetIndices[i][j] = randomRow;
                } else {
                    int randomCol = rand.nextInt(numCols);
                    subsetIndices[i][j] = randomCol;
                }
            }
        }
        return subsetIndices;
    }







//    /**
//     * Calculates the mean squared residue score for a given row of the data array
//     * @param the row for which to calculate the mean squared residue score
//     * @return the mean squared residue score for the given row
//     */


    double[] calculateMeanSquaredResidue(double[][] data, int[][] subsetIndices) {
        double[] meanSquaredResidueScores = new double[subsetIndices.length];
        for (int k = 0; k < subsetIndices.length; k++) {
            int[] rows = new int[6];
            int[] cols = new int[4];

            // Extract the row and column indices
            // Extract the row and column indices from the subsetIndices matrix
            for (int i = 0; i < 6; i++) {
                rows[i] = subsetIndices[k][i];
            }
            for (int i = 0; i < 4; i++) {
                cols[i] = subsetIndices[k][i + 6];
            }

            // Calculate the mean squared residue score for the submatrix
            double H = 0.0;
            int numEntries = rows.length * cols.length;
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < cols.length; j++) {
                    double value = data[rows[i]][cols[j]];
                    double sum = 0.0;
                    for (int ii = 0; ii < rows.length; ii++) {
                        for (int jj = 0; jj < cols.length; jj++) {
                            double diff = data[rows[ii]][cols[jj]] - value - data[rows[i]][cols[jj]] + data[rows[ii]][cols[j]];
                            sum += diff * diff;
                        }
                    }
                    H += sum;
                }
            }
            H = H / numEntries;

            // Store the mean squared residue score in the array
            meanSquaredResidueScores[k] = H;
        }

        return meanSquaredResidueScores;
    }



    private static double calculateMeanSquaredResidueForOneSubMatrix(double[][] data, double[] subsetValues) {
        int[] rows = new int[6];
        int[] cols = new int[4];

        // Extract the row and column indices from the subsetValues array
        for (int i = 0; i < 6; i++) {
            rows[i] = (int) subsetValues[i];
        }
        for (int i = 0; i < 4; i++) {
            cols[i] = (int) subsetValues[i + 6];
        }

        // Calculate the mean squared residue score for the submatrix
        double H = 0.0;
        int numEntries = rows.length * cols.length;
        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < cols.length; j++) {
                double value = data[rows[i]][cols[j]];
                double sum = 0.0;
                for (int ii = 0; ii < rows.length; ii++) {
                    for (int jj = 0; jj < cols.length; jj++) {
                        double diff = data[rows[ii]][cols[jj]] - value - data[rows[i]][cols[jj]] + data[rows[ii]][cols[j]];
                        sum += diff * diff;
                    }
                }
                H += sum;
            }
        }
        H = H / numEntries;

        return H;
    }

    public static double calculateMeanSquaredResidue(double[] alfa) {
        double H = 0.0;
        int numEntries = 24; // 6 rows x 4 columns

        // Calculate the mean squared residue score for the submatrix
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 4; j++) {
                double value = alfa[i];
                double sum = 0.0;
                for (int ii = 0; ii < 6; ii++) {
                    for (int jj = 0; jj < 4; jj++) {
                        double diff = alfa[ii] - value - alfa[i] + alfa[ii];
                        sum += diff * diff;
                    }
                }
                H += sum;
            }
        }
        H = H / numEntries;

        return H;
    }



    double calculateMeanSquaredResidueForOneSubMatrix(double[][] data, int[] subsetIndices) {
        // Extract the row and column indices from the subsetIndices array
        int[] rows = new int[6];
        int[] cols = new int[4];

        for (int i = 0; i < 6; i++) {
            rows[i] = subsetIndices[i];
        }
        for (int i = 0; i < 4; i++) {
            cols[i] = subsetIndices[i + 6];
        }

        // Calculate the mean squared residue score for the submatrix
        double H = 0.0;
        int numEntries = rows.length * cols.length;
        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < cols.length; j++) {
                double value = data[rows[i]][cols[j]];
                double sum = 0.0;
                for (int ii = 0; ii < rows.length; ii++) {
                    for (int jj = 0; jj < cols.length; jj++) {
                        double diff = data[rows[ii]][cols[jj]] - value - data[rows[i]][cols[jj]] + data[rows[ii]][cols[j]];
                        sum += diff * diff;
                    }
                }
                H += sum;
            }
        }
        H = H / numEntries;

        return H;
    }




    double[][] XX(int[][] subsetIndices, double[][] data) {
        int numRows = subsetIndices.length;
        int numCols = subsetIndices[0].length;
        double[][] subsetValues = new double[numRows][numCols];

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                int rowIndex = subsetIndices[i][j] / data[0].length;
                int colIndex = subsetIndices[i][j] % data[0].length;
                subsetValues[i][j] = data[rowIndex][colIndex];
            }
        }

        return subsetValues;
    }


    void init(){

        String filePath = "/home/skstar/Desktop/HuntingTon_500*50.xlsx";
        DataAndSize dataAndSize = mainData(filePath);

        double[][] data = dataAndSize.getData();
        int numRows = dataAndSize.getNumRows();
        int numCols = dataAndSize.getNumCols();



        int[][] subsetIndices = selectSubsetIndices(numRows,numCols);
        double[][] XX = XX(subsetIndices,data);



        XX=sort_and_index(XX);
        for(int i=0;i<D;i++) {
            alfa[i]=XX[0][i];
        }
        for(int i=0;i<D;i++) {
            beta[i]=XX[1][i];
        }
        for(int i=0;i<D;i++) {
            delta[i]=XX[2][i];
        }

    }










    double[][] solution() {

        init();


        int iter=1;
        double previousBestFitness = Double.MAX_VALUE;
        while(iter<=maxiter) {

            for(int j=0;j<D;j++) {
                a[j]=2.0-((double)iter*(2.0/(double)maxiter));
            }

            for(int i=0;i<N;i++)
            {
                for(int j=0;j<D;j++)
                {
                    r1=Math.random();
                    r2=Math.random();
                    for(int ii=0;ii<D;ii++) {
                        A1[ii]=2.0*a[ii]*r1-a[ii];
                    }
                    for(int ii=0;ii<D;ii++) {
                        C1[ii]=2.0*r2;
                    }

                    X1[i][j]=alfa[j]-A1[j]*(Math.abs(C1[j]*alfa[j]-XX[i][j]));
//                    X1= specificBounds(X1);

                    r1=Math.random();
                    r2=Math.random();
                    for(int ii=0;ii<D;ii++) {
                        A2[ii]=2.0*a[ii]*r1-a[ii];
                    }
                    for(int ii=0;ii<D;ii++) {
                        C2[ii]=2.0*r2;
                    }

                    X2[i][j]=beta[j]-A2[j]*(Math.abs(C2[j]*beta[j]-XX[i][j]));
//                    X2= specificBounds(X2);

                    r1=Math.random();
                    r2=Math.random();
                    for(int ii=0;ii<D;ii++) {
                        A3[ii]=2.0*a[ii]*r1-a[ii];
                    }
                    for(int ii=0;ii<D;ii++) {
                        C3[ii]=2.0*r2;
                    }

                    X3[i][j]=delta[j]-A3[j]*(Math.abs(C3[j]*delta[j]-XX[i][j]));
//                    X3= specificBounds(X3);


                    XX[i][j]=(X1[i][j]+X2[i][j]+X3[i][j])/3.0;

                }
            }
//            XX= specificBounds(XX);
            XX=sort_and_index(XX);

            for(int i=0;i<D;i++) {
                XX[N-1][i]=XX[0][i];
            }

            for(int i=0;i<D;i++) {
                alfa[i]=XX[0][i];
            }
            for(int i=0;i<D;i++) {
                beta[i]=XX[1][i];
            }
            for(int i=0;i<D;i++) {
                delta[i]=XX[2][i];
            }

//            for (int i = 0; i < alfa.length; i++) {
//                double value = alfa[i];
//                System.out.println("Value at index " + i + ": " + value);
//            }


            double bestFitness = calculateMeanSquaredResidue(alfa);
//                BESTVAL = new double[][] { { bestFitness }, alfa };
//            System.out.println(iter+" no . "+bestFitness);
            if (bestFitness < previousBestFitness) {
                System.out.println(bestFitness);
                previousBestFitness = bestFitness;
            } else {
                System.out.println(previousBestFitness);
            }

//            System.out.println(bestFitness);

            iter++;
        }




        double[][] bestSoliutionFound=new double[2][D];
        for(int i=0;i<D;i++) {
            bestSoliutionFound[1][i]=alfa[i];
        }
        bestSoliutionFound[0][0]=calculateMeanSquaredResidue(alfa);
        return bestSoliutionFound;

    }








    public double getOptimizeSoliution() {
        double[][] in=solution();
        return in[0][0];
//        System.out.println("Optimized value = "+in[0][0]);
//        for(int i=0;i<D;i++) {
//            System.out.println("x["+i+"] = "+in[1][i]);
//        }
    }

}


