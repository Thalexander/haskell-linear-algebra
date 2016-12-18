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

nearEqual a b = a - b >= -10**(-10) && a-b <= 10**(-10)

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

type LinearSystem = [Plane]

swap' :: LinearSystem -> Plane -> Plane -> LinearSystem
swap' [] _ _ = []
swap' (x:xs) a b
    | x == a = b : swap' xs a b
    | x == b = a : swap' xs a b
    | otherwise = x : swap' xs a b

swap :: LinearSystem -> Int -> Int -> LinearSystem
swap m a b = swap' m (m !! a) (m !! b)

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
