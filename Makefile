CPU_CC        := $(CC)
CPU_OPTS      := --std=c89 -Wall -Wno-format -Wno-unused-function -O3 # -Werror
CPU_MODULES   := exe lib bundled/klib bundled/xxhash bundled/xoshiro
CPU_LIBS	  := -lm -lz
CPU_BUILD_DIR := build
CPU_OBJ_DIR   := $(CPU_BUILD_DIR)/objs
CPU_SRC_DIR   := $(foreach item,$(CPU_MODULES),$(item)/src)
CPU_INC_DIR   := $(foreach item,$(CPU_MODULES),$(item)/include)
CPU_SRC       := $(foreach sdir,$(CPU_SRC_DIR),$(wildcard $(sdir)/*.c))
CPU_INC       := $(addprefix -I,$(CPU_INC_DIR))
CPU_BASE      := $(foreach sdir,$(CPU_SRC_DIR),$(notdir $(wildcard $(sdir)/*.c)))
CPU_OBJ       := $(patsubst %.c,$(CPU_OBJ_DIR)/%.o,$(CPU_BASE))
# SANITIZER    := -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer

.PHONY: all checkdirs clean

all: checkdirs build/solidclust

build/solidclust: $(CPU_OBJ)
	$(CPU_CC) $(SANITIZER) $^ -o $@ $(CPU_LIBS)

checkdirs: $(CPU_OBJ_DIR)

$(CPU_OBJ_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(CPU_BUILD_DIR) 

.SECONDEXPANSION:
PERCENT = %

$(CPU_OBJ): %.o : $$(filter $$(PERCENT)/$$(notdir %).c, $(CPU_SRC))
	$(CPU_CC) $(CPU_OPTS) $(SANITIZER) $(CPU_INC) -c $< -o $@
