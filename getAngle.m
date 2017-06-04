function theta = getAngle(oppositeLength, adjacentLength1, adjacentLength2)
theta = acos((adjacentLength1^2 + adjacentLength2^2 - oppositeLength^2)/(2*adjacentLength1*adjacentLength2));