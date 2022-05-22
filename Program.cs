using Linal;

class Program
{
    static public void Main(String[] args)
    {
        // Простая система 3х3
        double[,] a = {
            {5, 3, -2},
            {2, 1, -1},
            {3, -2, -3},
        };
        var A = new Matrix(a);
        var v1 = new Vector(-1, 0, 2);

        Console.WriteLine(A.ToString() + "----\n" + v1.ToString());
        Console.WriteLine("Roots (GaussJordan)   : " + Solver.GaussJordan(A, v1).ToString());
        Console.WriteLine("Roots (MatrixSolution): " + Solver.MatrixSolution(A, v1).ToString());
        Console.WriteLine();

        // Несовместная система
        double[,] b = {
            {2, 3, -1, 1},
            {8, 12, -9, 8},
            {4, 6, 3, -2},
            {2, 3, 9, -7},
        };
        var B = new Matrix(b);
        var v2 = new Vector(1, 3, 3, 3);

        Console.WriteLine(B.ToString() + "----\n" + v2.ToString());
        Console.WriteLine("Roots (GaussJordan)   : " + Solver.GaussJordan(B, v2).ToString());
        Console.WriteLine("Roots (MatrixSolution): " + Solver.MatrixSolution(B, v2).ToString());
        Console.WriteLine();

        // Совместная система
        double[,] c = {
            {2, 5, 4, 1},
            {1, 3, 2, 1},
            {2, 10, 9, 7},
            {3, 8, 9, 2},
        };
        var C = new Matrix(c);
        var v3 = new Vector(20, 11, 40, 37);

        Console.WriteLine(C.ToString() + "----\n" + v3.ToString());
        Console.WriteLine("Roots (GaussJordan)   : " + Solver.GaussJordan(C, v3).ToString());
        Console.WriteLine("Roots (MatrixSolution): " + Solver.MatrixSolution(C, v3).ToString());
        Console.WriteLine();
    }
}
