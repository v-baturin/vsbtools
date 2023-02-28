function matrix = VectoMat(vec)

matrix(1,1) = vec(1);
matrix(2,2) = vec(2);
matrix(3,3) = vec(3);
matrix(2,1) = vec(4);
matrix(1,2) = 0;
matrix(3,1) = vec(5);
matrix(1,3) = 0;
matrix(3,2) = vec(6);
matrix(2,3) = 0;
 
