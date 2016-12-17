module LinearAlgebra

where

newtype Vector = Vector (Double,Double,Double)
    deriving (Eq,Ord, Read, Show)
    
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
orthogonal a b = (dot a b) == 0
