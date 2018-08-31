#!/bin/bash

number=0
for i in `seq 1 49` ; do
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	let number=$number+1 ;./a.out $number &
	sleep 60s
done
