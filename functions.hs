module LinearAlgebra 

where

newtype Vector = Vector (Double,Double,Double)
    deriving (Eq, Ord, Read, Show)
    
instance Num Vector where
    (Vector (a,b,c)) + (Vector (x,y,z)) = Vector (a+x,b+y,c+z) 
    (Vector (a,b,c)) - (Vector (x,y,z)) = Vector (a-x,b-y,c-z)
    (Vector (a,b,c)) * (Vector (x,y,z)) = Vector (b*z - c*y, -(a*z - c*x), a*y - b*x)
    negate (Vector (a,b,c)) = Vector (-a,-b,-c)
    abs (Vector (a,b,c)) = Vector (abs a, abs b, abs c)
    signum (Vector (a,b,c)) = Vector (signum a, signum b, signum c)
    fromInteger a = Vector (fromInteger a, fromInteger a, fromInteger a)
    
scalarMult :: Double -> Vector -> Vector
scalarMult k (Vector (a,b,c)) = Vector (k*a, k*b, k*c)

magnitude :: Vector -> Double
magnitude (Vector (a, b, c)) = (a ** 2 + b ** 2 + c ** 2) ** 0.5

normalise :: Vector -> Vector
normalise v@(Vector (a,b,c))
    | mag == 0 = error "Cannot normalise the zero vector"
    | otherwise = Vector (a/mag, b/mag, c/mag)
    where mag = magnitude v
    
dot :: Vector -> Vector -> Double
dot (Vector (a,b,c)) (Vector (x,y,z)) = a*x + b*y + c*z

angle :: Vector -> Vector -> Double
angle a b = acos $ (dot a b) / ((magnitude a) * (magnitude b))

parallel :: Vector -> Vector -> Bool
parallel a b
    | magnitude a == 0 && magnitude b == 0 = True
    | magnitude a == 0 || magnitude b == 0 = False
    | otherwise = normalise a == normalise b
    
orthogonal :: Vector -> Vector -> Bool
orthogonal a b = nearEqual (dot a b) 0

componentParallel :: Vector -> Vector -> Vector
componentParallel v basis = scalarMult (dot v (normalise basis)) (normalise basis)

componentOrthogonal :: Vector -> Vector -> Vector
componentOrthogonal v basis = (v -) $ componentParallel v basis

areaParallelogram :: Vector -> Vector -> Double
areaParallelogram a b = magnitude $ a * b

areaTriangle :: Vector -> Vector -> Double
areaTriangle a b = (0.5 *) . magnitude $ a * b

first :: Vector -> Double
first (Vector (a,_,_)) = a

-- An alternative to == which accounts for floating point errors
nearEqual :: Double -> Double -> Bool
nearEqual a b = a - b >= -10**(-10) && a-b <= 10**(-10)

notNearZero :: Double -> Bool
notNearZero a = not $ nearEqual a 0

data Plane = Plane Vector Double
    deriving (Ord, Read, Show)
    
instance Eq Plane where
    p@(Plane f m) == q@(Plane g n) 
        | f == (Vector (0,0,0)) && g == (Vector (0,0,0)) = m == n
        | f == (Vector (0,0,0)) || g == (Vector (0,0,0)) = False
        | otherwise = orthogonal v f && orthogonal v g
            where Vector (a,b,c) = getPoint p
                  Vector (x,y,z) = getPoint q
                  v       = Vector (a-x,b-y,c-z)

on :: Vector -> Plane -> Bool
on v (Plane direction k) = (dot v direction) == k

firstNonZero :: Vector -> Maybe (Double, Int)
firstNonZero (Vector (a,b,c))
    | a /= 0 = Just (a, 0)
    | b /= 0 = Just (b, 1)
    | c /= 0 = Just (c, 2)
    | otherwise = Nothing
    
getPoint :: Plane -> Vector
getPoint (Plane v d)
        | firstNonZero v == Nothing = error "No valid point"
        | snd nonZero == 0 = Vector (d/coefficient,0,0)
        | snd nonZero == 1 = Vector (0,d/coefficient,0)
        | snd nonZero == 2 = Vector (0,0,d/coefficient)
    where Just nonZero = firstNonZero v
          coefficient = fst nonZero
          
planeParallel :: Plane -> Plane -> Bool
planeParallel (Plane v _) (Plane q _) = parallel v q

multPlane :: Plane -> Double -> Plane
multPlane (Plane v d) mult = Plane (scalarMult mult v) (d * mult)

addPlanes :: Plane -> Plane -> Plane
addPlanes (Plane v d) (Plane q e) = Plane (v + q) (d+ e)

-- A list of Planes that constructs a system so it can be used to find the solutions to where the planes meet
type LinearSystem = [Plane]

-- Helper function for swap
swap' :: LinearSystem -> Plane -> Plane -> LinearSystem
swap' [] _ _ = []
swap' (x:xs) a b
    | x == a = b : swap' xs a b
    | x == b = a : swap' xs a b
    | otherwise = x : swap' xs a b

-- Swap two rows in the LinearSystem
swap :: LinearSystem -> Int -> Int -> LinearSystem
swap m a b = swap' m (m !! a) (m !! b)

-- Multiply a row in the LinearSystem
multRow :: LinearSystem -> Double -> Int -> LinearSystem
multRow system multiplier row = map (\x -> if x == rowPlane then multPlane x multiplier else x) system
    where rowPlane = system !! row

addRows :: LinearSystem -> Int -> Int -> LinearSystem
addRows system a b = map (\x -> if x == rowA then addPlanes rowA rowB else x) system
    where rowA = system !! a
          rowB = system !! b
          
addPlaneToRow :: LinearSystem -> Plane -> Int -> LinearSystem
addPlaneToRow system plane row = map (\x -> if x == rowPlane then addPlanes x plane else x) system
    where rowPlane = system !! row

addMultiple :: LinearSystem -> Double -> Int -> Int -> LinearSystem
addMultiple system mult a b = addPlaneToRow system plane a
    where plane = multPlane (system !! b) mult
    
firstNonZeroRow :: LinearSystem -> Int -> Int -> Maybe Int
firstNonZeroRow [] _ _ = Nothing
firstNonZeroRow system@((Plane v d):xs) row pos
    | nonZero == Nothing = Nothing
    | (snd value) /= pos = firstNonZeroRow xs (row + 1) pos
    | otherwise = Just row
    where nonZero    = firstNonZero v
          Just value = nonZero

-- Swaps rows around to ensure the pos in the LinearSystem is non zero
makeNonZeroRow :: LinearSystem -> Int -> Int -> Maybe LinearSystem
makeNonZeroRow system row pos
    | firstRow == Nothing = Nothing
    | otherwise = Just $ swap system row swapRow
    where firstRow = firstNonZeroRow system row pos
          Just swapRow = firstRow

getVectorPos :: Vector -> Int -> Double
getVectorPos (Vector (a,b,c)) pos
    | pos == 0 = a
    | pos == 1 = b
    | pos == 2 = c
    | otherwise = error "Invalid position"

coefficient :: LinearSystem -> Int -> Int -> Double
coefficient system row pos = getVectorPos v pos
    where (Plane v _) = system !! row 

-- Reduce the values below the first row in a position to 0, returns the LinearSystem without the first row
clearBelow :: LinearSystem -> Int -> LinearSystem
clearBelow (_:[]) _ = [] 
clearBelow system@(x:(y:ys)) pos = (addMultiple system mult 1 0 !! 1) : clearBelow (x:ys) pos
    where coX  = coefficient system 0 pos
          coY  = coefficient system 1 pos
          mult = - (coY / coX)
    
-- Alternate version of clearBelow that returns with the first row
reduceRowsBelow :: LinearSystem -> Int -> LinearSystem
reduceRowsBelow system@(x:_) pos = x : clearBelow system pos

reduceColumn :: LinearSystem -> Int -> LinearSystem
reduceColumn [] _ = []
reduceColumn system@(x:_) pos
    | pos >= 3 = system
    | nonZeroRow == Nothing = reduceColumn system (pos + 1)
    | nearEqual c 0 = reduceColumn newSystem pos
    | otherwise = x : reduceColumn (clearBelow system pos) (pos + 1)
    where c = coefficient system 0 pos
          nonZeroRow     = makeNonZeroRow system 0 pos
          Just newSystem = nonZeroRow
    
-- Converts LinearSystem to triangular form
triangularForm :: LinearSystem -> LinearSystem
triangularForm system = reduceColumn system 0

nonZeroValues :: Plane -> Bool
nonZeroValues (Plane v _)
    | firstNonZero v == Nothing = False
    | otherwise                 = True

-- Change leading values in each row to 1 and reduces values in the same position below it to 0
simplifySystem :: LinearSystem -> LinearSystem
simplifySystem []     = []
simplifySystem system@(x:[]) = (multPlane x (1 / coX)):[]
    where Plane v _   = x
          Just index' = firstNonZero v
          index       = snd index'
          coX         = coefficient system 0 index
simplifySystem system@(x:xs)
    | nonZeroValues x == False = x : simplifySystem xs
    | otherwise                = (multPlane x (1 / coX)) : (simplifySystem $ clearBelow system index)
    where Plane v _   = x
          Just index' = firstNonZero v
          index       = snd index'
          coX         = coefficient system 0 index

-- Converts to reduced row echelon format
rref :: LinearSystem -> LinearSystem
rref system = reverse $ simplifySystem $ triangular
    where triangular = reverse $ triangularForm system
 
-- Number of solutions to a LinearSystem, if Unique has the vector that satisfies it, in future the Infinite solution could contain the parameterised solution   
data SystemSolution = Infinite | Unique Vector | None
    deriving (Eq, Ord, Show, Read)

-- Finds if a LinearSystem in RREF has contradictory equations of the form 0 = k
contradictoryEquations :: LinearSystem -> Bool
contradictoryEquations [] = False
contradictoryEquations (x@(Plane _ d):xs)
    | nonZeroValues x == False && not (nearEqual d 0) = True
    | otherwise                                       = contradictoryEquations xs
    
pivots :: LinearSystem -> Int
pivots [] = 0
pivots (x@(Plane v _):xs)
    | firstNonZero v == Nothing = 0 + pivots xs
    | otherwise                 = 1 + pivots xs

-- Finds if a LinearSystem in RREF has too few pivots to have a unique solution
tooFewPivots :: LinearSystem -> Bool
tooFewPivots = (< 3) . pivots

constantTerm :: Plane -> Double
constantTerm (Plane _ d) = d

-- Helper function for solveSystem, requires input to be in RREF
findUniqueSolution :: LinearSystem -> Vector
findUniqueSolution (x:(y:(z:_))) = Vector (constantTerm x, constantTerm y, constantTerm z)
 
-- Finds if there are infinite solutions, none, or a unique solution to a LinearSystem   
solveSystem :: LinearSystem -> SystemSolution
solveSystem system 
    | contradictoryEquations rrefSystem = None
    | tooFewPivots rrefSystem           = Infinite
    | otherwise                         = Unique $ findUniqueSolution rrefSystem
    where rrefSystem = rref system
