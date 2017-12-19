public class MatrixCalculationImpl implements MatrixCalculation {

    private int formula;
    private int orderOfFunction;
    private double C;
    private double minError;

    MatrixCalculationImpl(int formula, int order, double c) {
        this.formula = formula;
        this.orderOfFunction = order;
        this.C = c;
        this.minError = 0.001;
    }

    public double[][] calculateMatrix(double betta, double gamma) {
        return null;
    }

    @Override
    public double[][] transposePMatrix(double[][] matrix) {
        return null;
    }

    @Override
    public double[][] multiplyPMatrixByPTransposed(double[][] matrix) {
        return null;
    }
}

