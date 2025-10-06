<template>
  <button v-show="showBackToTop" class="back-to-top" @click="handleClick" aria-label="回到顶部"
    @touchstart="handleTouchStart" @touchmove="handleTouchMove" @touchend="handleTouchEnd"
    :style="{ top: buttonTop + 'px' }">
    <i class="fas fa-arrow-up"></i>
  </button>
</template>

<script>
export default {
  name: 'BackToTop',
  data() {
    return {
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
    window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: this.buttonTop } }))
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.handleScroll)
    window.removeEventListener('floating-buttons-base-top', this.syncBaseTop)
  },
  methods: {
    // 接收“下方按钮”的基准位置（BackToTop 自己就是下方按钮），直接对齐
    syncBaseTop(e) {
      const base = e?.detail?.baseTop
      if (typeof base === 'number') this.buttonTop = base
    },
    handleScroll() {
      if (!this.isDragging) {
        this.showBackToTop = window.scrollY > 180
        // 当置顶按钮可见时，持续派发基准位置，保证两按钮间距始终一致
        if (this.showBackToTop) {
          window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: this.buttonTop } }))
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
      this.isDragging = true
      this.touchMoved = false
      this.startY = e.touches[0].clientY
      this.initialTop = this.buttonTop
    },
    handleTouchMove(e) {
      this.touchMoved = true

      const currentY = e.touches[0].clientY
      const diffY = currentY - this.startY

      let newTop = this.initialTop + diffY

      const buttonHeight = 40
      const GAP = 48
      // 下方按钮需要预留 GAP 的上方空间，避免与 TOC 重合
      const minTop = 20 + GAP
      const maxTop = window.innerHeight - buttonHeight - 20

      newTop = Math.max(minTop, newTop)
      newTop = Math.min(maxTop, newTop)

      this.buttonTop = newTop

      // BackToTop 为下方按钮：派发自身作为基准位置
      window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: newTop } }))

      e.preventDefault()
    },
    handleTouchEnd() {
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
  background-color: #fff; /* 纯白 */
  color: var(--bs-body-color);
  border: 1px solid var(--bs-border-color);
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}

.back-to-top:hover {
  background-color: #fff; /* hover 仍保持纯白 */
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.18);
}

.back-to-top:active {
  transform: scale(0.95);
}
</style>