public interface MatrixCalculation {
    double[][] calculateMatrix(double betta, double gamma);
    double[][] transposePMatrix(double[][] matrix);
    double[][] multiplyPMatrixByPTransposed(double[][] matrix);
}
