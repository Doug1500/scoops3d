package main

import (
    "bufio"
    "fmt"
    "log"
    "os"
)

func main() {
    file, err := os.Open("A_emb20_out.txt")
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()

    important:= "1st"
    important2:= "2nd"
    counter := 0.0

    scanner := bufio.NewScanner(file)
    for scanner.Scan() {
        a:=scanner.Text()
        if counter < 2.0{
            if a == "3D POTENTIAL FAILURE"{
                important = a
            }
            if important == "3D POTENTIAL FAILURE"{
                important2 = a
                counter = counter + 1.0
            }
        }
    }
    fmt.Println(important2)
    important3 := float64(important2)

    if err := scanner.Err(); err != nil {
        log.Fatal(err)
    }
}