<template>
  <button v-show="showBackToTop" class="back-to-top" @click="handleClick" aria-label="回到顶部"
    @touchstart.prevent.stop="handleTouchStart" @touchmove.prevent.stop="handleTouchMove" @touchend.prevent.stop="handleTouchEnd"
    :style="{ top: buttonTop + 'px' }">
    <i class="fas fa-arrow-up"></i>
  </button>
</template>

<script>
export default {
  name: 'BackToTop',
  data() {
    return {
      // 事件源标识与 raf 派发节流
      sourceId: 'btt',
      rafPending: false,
      rafLastBaseTop: null,

      showBackToTop: false,
      isDragging: false,
      startY: 0,
      initialTop: 0,
      buttonTop: window.innerHeight - 100,
      touchMoved: false
    }
  },
  mounted() {
    window.addEventListener('scroll', this.handleScroll)
    window.addEventListener('floating-buttons-base-top', this.syncBaseTop)
    this.handleScroll()
    this.buttonTop = window.innerHeight - 100
    // 初始同步：派发当前作为“下方按钮”的基准位置
    this.rafDispatchBaseTop(this.buttonTop)
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.handleScroll)
    window.removeEventListener('floating-buttons-base-top', this.syncBaseTop)
  },
  methods: {
    // 工具：统一边界
    getBounds() {
      const BUTTON_HEIGHT = 40
      const GAP = 48
      const MARGIN = 20
      return {
        gap: GAP,
        minTop: MARGIN + GAP, // 作为下方按钮，需要为上方按钮预留 GAP
        maxTop: window.innerHeight - BUTTON_HEIGHT - MARGIN
      }
    },
    clampTop(top) {
      const { minTop, maxTop } = this.getBounds()
      return Math.max(minTop, Math.min(maxTop, top))
    },
    // raf 节流派发，避免频繁互相抖动
    rafDispatchBaseTop(baseTop) {
      this.rafLastBaseTop = baseTop
      if (this.rafPending) return
      this.rafPending = true
      requestAnimationFrame(() => {
        window.dispatchEvent(new CustomEvent('floating-buttons-base-top', {
          detail: { baseTop: this.rafLastBaseTop, source: this.sourceId }
        }))
        this.rafPending = false
      })
    },
    // 接收“下方按钮”的基准位置（BackToTop 自己就是下方按钮），对齐但不回收自身派发
    syncBaseTop(e) {
      const base = e?.detail?.baseTop
      const source = e?.detail?.source
      if (source === this.sourceId) return
      if (typeof base === 'number') this.buttonTop = this.clampTop(base)
    },
    handleScroll() {
      if (!this.isDragging) {
        this.showBackToTop = window.scrollY > 180
        // 当置顶按钮可见时，通过 raf 节流派发基准位置
        if (this.showBackToTop) {
          this.rafDispatchBaseTop(this.buttonTop)
        }
      }
    },
    backToTop() {
      window.scrollTo({ top: 0, behavior: 'smooth' })
    },
    handleClick() {
      if (!this.isDragging) {
        this.backToTop()
      }
    },
    handleTouchStart(e) {
      e.preventDefault()
      this.isDragging = true
      this.touchMoved = false
      this.startY = e.touches[0].clientY
      this.initialTop = this.buttonTop
    },
    handleTouchMove(e) {
      this.touchMoved = true

      const currentY = e.touches[0].clientY
      const diffY = currentY - this.startY

      let newTop = this.clampTop(this.initialTop + diffY)

      this.buttonTop = newTop

      // BackToTop 为下方按钮：派发自身作为基准位置
      this.rafDispatchBaseTop(newTop)

      e.preventDefault()
    },
    handleTouchEnd(e) {
      e.preventDefault()
      this.isDragging = false

      if (!this.touchMoved) {
        this.backToTop()
      }
    }
  }
}
</script>

<style scoped>
.back-to-top {
  position: fixed;
  right: 30px;
  width: 40px;
  height: 40px;
  background-color: var(--app-btn-bg);
  color: var(--app-text);
  border: 1px solid var(--app-border);
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: var(--app-btn-shadow);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}

.back-to-top:hover {
  background-color: var(--app-btn-hover-bg);
  box-shadow: var(--app-btn-hover-shadow);
}

.back-to-top:active {
  transform: scale(0.95);
}
</style>