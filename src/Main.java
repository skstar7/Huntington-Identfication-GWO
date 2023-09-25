import java.util.Arrays;

public class Main {
    public static double mean(double arr[]){
        double total = 0;
        for (int i = 0 ;i<arr.length;i++){
            total +=arr[i];
        }

        return total/arr.length;
    }


    // Function for calculating median
    public static double findMedian(double arr[]) {
        // First we sort the array
        Arrays.sort(arr);
        int n = arr.length;
        // check for even case
        if (n % 2 != 0)
            return arr[n / 2];

        return (arr[(n - 1) / 2] + arr[n / 2]) / 2.0;
    }

    //method to calculate standard deviation
    public static double standardDeviation(double arr[]) {

        // calculate the standard deviation
        double standardDeviation = 0.0;
        for (double num : arr) {
            standardDeviation += Math.pow(num - mean(arr), 2);
        }

        return Math.sqrt(standardDeviation /arr.length);
    }
    public static double minCalculate(double arr[]){
        double min = arr[0];
        //Loop through the array
        for (int i = 0; i < arr.length; i++) {
            //Compare elements of array with min
            if(arr[i] <min)
                min = arr[i];
        }
        return min;

    }
    public static double maxCalculate(double arr[]){
        double max  = arr[0];
        for (int i = 0; i < arr.length; i++) {
            //Compare elements of array with min
            if(arr[i] >max)
                max= arr[i];
        }
        return  max;
    }
    public static void main(String[] args) {
        int maxIterations = 200; // Set the desired maximum number of iterations

        double[] soln = new double[50];
        grey_Wolf_Optimization obj1 = new grey_Wolf_Optimization(maxIterations);
        for (int i = 0; i < 50; i++) {

            System.out.println("for " + i + " no. run ");
            soln[i] = obj1.getOptimizeSoliution();
            System.out.println(soln[i]);
        }
        System.out.println("Mean of the sphear function is " + mean(soln));
        System.out.println("Median of the sphear function is " + findMedian(soln));
        System.out.println("Standard deviation of the sphear function is "+standardDeviation(soln));
        System.out.println("Min of the sphere function " + minCalculate(soln));
        System.out.println("Max of the sphere function " + maxCalculate(soln));




        // You can also access the optimized solution's variables if needed
//         double[] optimizedSolution = optimizer.getOptimizedSolution();
//         System.out.println("Optimized solution: " + Arrays.toString(optimizedSolution));
 }
}
