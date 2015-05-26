# General

A simple hangman-like game! Given a maximum number of guesses, the number of letters to remove from the word, and a word.
Play a short game of hangman with a random set of letters removed from the word.
Example game play provided below.

Note: Removing the same letter twice does not count as a removal, for instance hello as he__o is only removing 1 letter.

# Usage

```
mvn package
java -cp target/interview-hangman-0.0.1-SNAPSHOT.jar com.bina.interview.hangman.Hangman <max guess count> <missing letter count> <word>
```

# A winning game

```
java -jar target/interview-hangman.jar 8 3 helloworld

0 3 _e__owo__d
> k
:(

1 3 _e__owo__d
> l
:)

1 2 _ellowo_ld
> t
:(

2 2 _ellowo_ld
> t
:(

3 2 _ellowo_ld
> h
:)

3 1 hellowo_ld
> r
:)

helloworld
winner!
```

# A losing game

```
java -jar target/interview-hangman.jar 2 1 cat

0 1 c_t
> k
:(

0 1 c_t
> help
> j
:(

cat
loser!
```
