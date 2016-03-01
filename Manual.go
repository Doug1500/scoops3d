package main

import (
	"fmt"
	"scoops3d/Part"
)

func main() {

	fmt.Println("Welcome to GASlope 3D LEM")
	fmt.Println("Can you please x, y, z, R, alpha ")
	// FOS := 4.102000e+00   
	// Igood := []float64{561673, 6.11237e+06, 2835.46, 3179.28, 182.659}
	var x, y, z, R, alpha, th1, th2 float64
	fmt.Scan(&x, &y, &z, &R, &alpha)
	
	score, err := Part.Partmain(x, y, z, R, alpha*(180*7)/22)
	// score, err := Part.Partmain(x, y, z, R, alpha, FOS, Igood)
	fmt.Println(score, err, x, y, z, R, alpha*(180*7)/22 )

}