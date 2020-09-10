## Performance

[x] Pakiet BenchmarkTools posiada makro ``@btime`` , które wielokrotnie uruchamia kod i wyswietla rowniez alokacje pamieci,
``@benchmark`` wyswietla bardziej szczegolowe informacje

[x] Julia najpierw iteruje po wartosciach w kolumnach a dopiero potem po wierszach (tak umieszcza je w pamieci).

[x] ``@code_warntype baz()`` - wyswietla typy wystepujace w przebiegu funkcji baz()

[x]
```Julia
julia> @macroexpand x = @. A + b*c
:(x = (+).(A, (*).(b, c)))
```

[x] deklarowanie typow w funkcji nie jest koniecznie potrzebne kiedy dobrze jest napisany przebieg funckji, tak zeby komplilator wiedzial jakich typow uzyc. Jezeli ``@code_warntype`` zwroci unie typow przy wywolaniu oznacza to ze jest cos nie tak.

[x] Profiler - wskazuje dokladny przebieg funkcji z informacja jakie linijki byly najdluzej wykonywane

[x] Globalne wartosci deklarowac jako const (np. macierze ktore beda przekazywane do funkcji)

[x] Funkcje inlinowe ``@inline``

[x] Generowana funkcja na poziomie kompilacji kodu (74)

[x] Subnormal numbers - numery bardzo bliskie zeru zamieniane na zero

[x] Kiedy macierz jest np. transponowana, to typ takiej macierzy moze byc ``Adjoint`` czyli odpowiednik macierzy wejsciowej, to taka macierz najlepiej iterowac po wierszach

[x] Do inicjalizacji macierzy wykorzystywac funkcje ``fill()``

[x] ``@inbounds`` - wynikowy kod nie sprawdza czy odwolanie do tablicy nie zostalo przekroczone (a przeciez w wielu miejsach kodu wiemy ze sie nie przekroczy)

[x] Prealokacja macierzy oraz przekazywanie macierzy przez referencje do funkcji

[x] Wykorzysywac wersje funkcji z ``!`` ktora dziala na oryginalnej wersji tablicy, bez tworzenia kopii

[x] Do operacji grupowych na tablicach wykorzysywac wersje funkcji z kropka ``.``, np. ``A .+ B``, ``sin.(A)``

[x] ``@view`` pozwala podzielic macierz na czesci bez tworzenia kopii, np. ``@view(x[:,i])``

[x] ``@simd`` - wykorzystanie mozliwosci wykonywania operacji o podobnej strutkruze (waruneki: to co sie dzieje miedzy iteracjami musi byc niezalezne i w srodku nie moze byc wywolyania funkcji) - boucd checking jest domyslnie wylaczone. Jest tez biblioteka **SIMD.jl**

[x] Iterowac po obiektach a nie po indeksach ``for a in A``. A jezeli wartosci indeksu jest wymagany do oblczen to skorzystac z funkcji ``eachindex()``, np. ``for i in eachindex(A)``

[] Wielowatkowosc ``export JULIA_NUP_THREADS = 4``. ``@threads``

## Metoda potegowa:
- im wektor startowy jest bliższy poszukiwanemu wektorowi własnemu, tym mniej iteracji trzeba będzie wykonać.
