package main

import (
	"fmt"
	"os/exec"
	// "github.com/kardianos/osext"
)

func main() {
	err:= exec.Command("/home/yewintun/mygo/src/scoop3d/scoops3d/scoops3d.exe").Run()
	fmt.Printf("HAHAH",err)
	// filename, _ := osext.excutable()
	// fmt.Println(filename)
}