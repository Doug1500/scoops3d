package main

import (
	"fmt"
	"scoops3d/Part"
	"math"
)

func main() {

	fmt.Println("Welcome to GASlope 3D LEM")
	fmt.Println("Can you please x, y, z, R, alpha, th1, RI")

	suksd, phiksd, sukmean, phikmean, gamma, su11, phi11 := Part.SoilP()

	RITolerance := 0.01

	var x, y, z, R, alpha, th1, RI float64
	fmt.Scan(&x, &y, &z, &R, &alpha, &th1, &RI)

	// phi11[1] = phikmean[1] + phiksd*RI*math.Cos(th1)
	// su11[1] = sukmean[1] + suksd*RI*math.Sin(th1)
	
	phi11[0] = phikmean[0] + phiksd[0]*RI*math.Cos(th1)
	su11[0] = sukmean[0] + suksd[0]*RI*math.Sin(th1)

	FOSS, err := Part.Partmain(x, y, z, R, alpha*180*7/22, 0, su11, phi11, gamma)
	
	score := 10000.0
	if math.Abs(FOSS-1.0) < RITolerance && err == false {
		score = RI
	} else {
		score = RI + 100.0*math.Abs(FOSS-1.0) + 5.0
	}

	fmt.Println(phi11, su11, sukmean)
	fmt.Println(FOSS)
	fmt.Println(score, err)
	fmt.Println("RI:", RI)

}