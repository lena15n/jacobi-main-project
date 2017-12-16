import java.util.Arrays;

public class Test {

    private static final double C = 2;

    public static void main(String[] args) {
        //0 10 4 2 1 0.05
        double gamma = 2;
        double betta = 1;
        long time = System.currentTimeMillis();


        int order = 4;//10;
        double minAllowedError = 0.001;

        double[][] matrixDerivate = new MatrixCalculation(0, order, C, minAllowedError).calculateMatrix(betta, gamma);
        double[][] matrix = new MatrixCalculation(1, order, C, minAllowedError).calculateMatrix(betta, gamma);
        double[][] matrixIntegral = new MatrixCalculation(2, order, C, minAllowedError).calculateMatrix(betta, gamma);
//        print("derivate", matrixDerivate);
//        print("p", matrix);
//        print("integral", matrixIntegral);
        time = System.currentTimeMillis() - time;
        System.err.println("Time: " + time);
    }

    public static void print(String name, double[][] matrix) {
        System.err.println(name);
        int K = 1;
        for (double[] m : matrix) {
            System.err.print((K++) + " ");
            System.err.println(Arrays.toString(m));
        }
        System.err.println("===");
    }
}