``@code_llvm func(m)`` - pokazuje implementacje w LLVM

Julia najpierw iteruje po wartosciach w kolumnach a dopiero potem po wierszach (tak umieszcza je w pamieci). Nalezy to wziac gdy faktycznie bedziemy iterowac wartosci w macierzy.

Pakiet BenchmarkTools posiada makro ``@btime`` , które wielokrotnie uruchamia kod i wyswietla rowniez alokacje pamieci

``@code_warntype baz()`` - wyswietla typy wystepujace w przebiegu funkcji baz()

```Julia
julia> @macroexpand x = @. A + b*c
:(x = (+).(A, (*).(b, c)))
```
Rozpisuje jak wyglada fragment kodu po wykorzystaniu makra
