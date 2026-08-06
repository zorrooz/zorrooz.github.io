/**
 * 浮动按钮（BackToTop / TocDrawer）共享的定位与协调逻辑。
 * 两按钮通过 window CustomEvent 'floating-buttons-base-top' 协调纵向位置：
 * - mode 'match'：直接占据广播的 base
 * - mode 'stack'：挂在 base 上方（top = base - gap）
 * 支持触摸拖拽（constrained 到边界）与单击释放（触发 onRelease）。
 */
import { ref } from 'vue'

export interface FloatingButtonOptions {
  sourceId: string
  defaultTop: number
  mode: 'match' | 'stack'
  buttonHeight?: number
  gap?: number
  margin?: number
  /** 触摸结束且未发生拖拽时调用（单击语义，如回顶/打开抽屉） */
  onRelease: () => void
}

export function useFloatingButton(opts: FloatingButtonOptions) {
  const { sourceId, mode, onRelease } = opts
  const buttonHeight = opts.buttonHeight ?? 40
  const gap = opts.gap ?? 48
  const margin = opts.margin ?? 20

  const rafPending = ref(false)
  const rafLastBaseTop = ref<number | null>(null)
  const isDragging = ref(false)
  const startY = ref(0)
  const initialTop = ref(0)
  const buttonTop = ref(opts.defaultTop)
  const touchMoved = ref(false)

  const getBounds = () => ({
    gap,
    minTop: mode === 'stack' ? margin : margin + gap,
    maxTop:
      mode === 'stack'
        ? Math.max(0, window.innerHeight - buttonHeight - margin - gap)
        : window.innerHeight - buttonHeight - margin,
  })

  const clampTop = (top: number) => {
    const { minTop, maxTop } = getBounds()
    return Math.max(minTop, Math.min(maxTop, top))
  }

  /** 当前自身 top → 应该广播出去的 base */
  const baseFromCurrent = (top: number) => (mode === 'stack' ? top + gap : top)
  /** 收到的 base → 自身应该落位的 top */
  const topFromBase = (base: number) => (mode === 'stack' ? base - gap : base)

  function dispatchBaseTop() {
    rafLastBaseTop.value = baseFromCurrent(buttonTop.value)
    if (rafPending.value) return
    rafPending.value = true
    requestAnimationFrame(() => {
      if (rafLastBaseTop.value !== null) {
        window.dispatchEvent(
          new CustomEvent('floating-buttons-base-top', {
            detail: { baseTop: rafLastBaseTop.value, source: sourceId },
          }),
        )
      }
      rafPending.value = false
    })
  }

  function handleBaseTopEvent(e: Event) {
    const detail = (e as CustomEvent).detail
    const base = detail?.baseTop
    const source = detail?.source
    if (source === sourceId) return
    if (typeof base === 'number') buttonTop.value = clampTop(topFromBase(base))
  }

  function onTouchStart(e: TouchEvent) {
    e.preventDefault()
    isDragging.value = true
    touchMoved.value = false
    startY.value = e.touches[0].clientY
    initialTop.value = buttonTop.value
  }

  function onTouchMove(e: TouchEvent) {
    touchMoved.value = true
    const diffY = e.touches[0].clientY - startY.value
    buttonTop.value = clampTop(initialTop.value + diffY)
    dispatchBaseTop()
    e.preventDefault()
  }

  function onTouchEnd(e: TouchEvent) {
    e.preventDefault()
    isDragging.value = false
    if (!touchMoved.value) onRelease()
  }

  function subscribe() {
    window.addEventListener('floating-buttons-base-top', handleBaseTopEvent)
  }

  function unsubscribe() {
    window.removeEventListener('floating-buttons-base-top', handleBaseTopEvent)
  }

  return {
    isDragging,
    buttonTop,
    clampTop,
    dispatchBaseTop,
    onTouchStart,
    onTouchMove,
    onTouchEnd,
    subscribe,
    unsubscribe,
  }
}
