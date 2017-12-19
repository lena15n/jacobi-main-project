public class Tool {

    private static final double C = 2;
    private static final double MIN_ALLOWED_ERROR = 0.001;

    public static void main(String[] args) {
        double gamma = 2;
        double betta = 1;
        int order = 10;

        // set same MIN_ALLOWED_ERROR!!
//        Result derivateTime = runDerivate(gamma, betta, order);
//        Result pTime        = runP(gamma, betta, order);
//        Result integralTime = runIntegral(gamma, betta, order);

//        Result transposeTime = runTranspose(gamma, betta, order);
        Result multiplyTime  = runMultiplyPMatrixByPTransposed(gamma, betta, order);

//        System.out.println(String.format("derivateTime: %d", derivateTime.getTimeMillisParallel()));
//        System.out.println(String.format("pTime: %d", pTime.getTimeMillisParallel()));
//        System.out.println(String.format("integralTime: %d", integralTime.getTimeMillisParallel()));
//        System.out.println(String.format("transposeTime: %d", transposeTime.getTimeMillisParallel()));
        System.out.println(String.format("multiplyTime: %d", multiplyTime.getTimeMillisParallel()));
    }


    private static Result runDerivate(double gamma, double betta, int order) {
        return runJacobiCalc(0, gamma, betta, order);
    }

    private static Result runP(double gamma, double betta, int order) {
        return runJacobiCalc(1, gamma, betta, order);
    }

    private static Result runIntegral(double gamma, double betta, int order) {
        return runJacobiCalc(2, gamma, betta, order);
    }

    private static Result runJacobiCalc(int formula, double gamma, double betta, int order) {
        Result result = new Result();

//        long time = System.currentTimeMillis();
//        new MatrixCalculationImpl(formula, order, C).calculateMatrix(betta, gamma);
//        result.setTimeMillisSeries(System.currentTimeMillis() - time);

        long time = System.currentTimeMillis();
        new MatrixParallelCalculationImpl(formula, order, C).calculateMatrix(betta, gamma);
        result.setTimeMillisParallel(System.currentTimeMillis() - time);

        return result;
    }

    private static Result runTranspose(double gamma, double betta, int order) {
        Result result = new Result();
        int p = 1;
//        long time = System.currentTimeMillis();
//        new MatrixCalculationImpl(p, order, C).transposePMatrix(betta, gamma);
//        result.setTimeMillisSeries(System.currentTimeMillis() - time);

        MatrixParallelCalculationImpl mc = new MatrixParallelCalculationImpl(p, order, C);
        double[][] matrix = mc.toRectangleMatrix(mc.calculateMatrix(betta, gamma));

        long time = System.currentTimeMillis();
        mc.transposePMatrix(matrix);
        result.setTimeMillisParallel(System.currentTimeMillis() - time);

        return result;
    }

    private static Result runMultiplyPMatrixByPTransposed(double gamma, double betta, int order) {
        Result result = new Result();
        int p = 1;
//        long time = System.currentTimeMillis();
//        new MatrixCalculationImpl(p, order, C).multiplyPMatrixByPTransposed(betta, gamma);
//        result.setTimeMillisSeries(System.currentTimeMillis() - time);


        MatrixParallelCalculationImpl mc = new MatrixParallelCalculationImpl(p, order, C);
        double[][] matrix = mc.toRectangleMatrix(mc.calculateMatrix(betta, gamma));

        long time = System.currentTimeMillis();
        mc.multiplyPMatrixByPTransposed(matrix);
        result.setTimeMillisParallel(System.currentTimeMillis() - time);

        return result;
    }
}
