package main

import (
	"fmt"
	// "math/rand"
	// "math"
	"scoops3d/Part"
)

func main() {

	x, y, z, R, alpha :=  561186.0, 5117810.0, 3590.33, 2001.1, 139.631

	_, _, _, _, gamma, su11, phi11 := Part.SoilP()

	su11[0], phi11[0] = 1000.0, 40.0

	score, err := Part.Partmain(x, y, z, R, alpha, su11, phi11, gamma)

	fmt.Println(score, err)
}


