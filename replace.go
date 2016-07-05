package main

import (
        "io/ioutil"
        "log"
        "strings"
        "fmt"
)

func main() {
        input, err := ioutil.ReadFile("emb20mat3D.txt")
        if err != nil {
                log.Fatalln(err)
        }

        lines := strings.Split(string(input), "\n")

        for i, _ := range lines {
                // if strings.Contains(line, "Scoops3D") {
                //         lines[i] = "LOL"
                // }
                fmt.Println(i, lines[i][0]) 
                
        }
        output := strings.Join(lines, "\n")
        err = ioutil.WriteFile("emb20mat3D.txt", []byte(output), 0644)
        if err != nil {
                log.Fatalln(err)
        }
}