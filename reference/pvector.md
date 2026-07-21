# A push/pop capable vector

An R6 class that implements a persistent vector with push and pop
operations.

## Public fields

- `values`:

  A vector of values

## Methods

### Public methods

- [`pvector$new()`](#method-pvector-new)

- [`pvector$print()`](#method-pvector-print)

- [`pvector$length()`](#method-pvector-length)

- [`pvector$pop()`](#method-pvector-pop)

- [`pvector$push()`](#method-pvector-push)

- [`pvector$clone()`](#method-pvector-clone)

------------------------------------------------------------------------

### Method `new()`

Create a pvector

#### Usage

    pvector$new(values = c())

#### Arguments

- `values`:

  A vector of values

#### Returns

A new pvector

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print pvector

#### Usage

    pvector$print()

#### Returns

NULL

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

Get length of pvector

#### Usage

    pvector$length()

#### Returns

An integer

------------------------------------------------------------------------

### Method `pop()`

Pop vector

#### Usage

    pvector$pop()

#### Returns

Popped value

------------------------------------------------------------------------

### Method `push()`

Push values

#### Usage

    pvector$push(pushed.values)

#### Arguments

- `pushed.values`:

  A vector of values

#### Returns

NULL

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    pvector$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
