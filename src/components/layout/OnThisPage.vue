<template>
  <nav class="on-this-page">
    <div class="otp-header">
      <span class="otp-title">本页目录</span>
    </div>

    <ul class="otp-list">
      <li v-for="(item, idx) in toc" :key="idx" :class="['otp-item', { active: activeId === item.id }]">
        <a class="otp-link" role="button" tabindex="0" @click.prevent="scrollToId(item.id)" @keydown.enter.prevent="scrollToId(item.id)">
          <span class="otp-text">{{ item.text }}</span>
        </a>
        <ul v-if="item.children && item.children.length" class="otp-sublist">
          <li v-for="(sub, sIdx) in item.children" :key="sIdx" :class="['otp-subitem', { active: activeId === sub.id }]">
            <a class="otp-sublink" role="button" tabindex="0" @click.prevent="scrollToId(sub.id)" @keydown.enter.prevent="scrollToId(sub.id)">
              <span class="otp-subtext">{{ sub.text }}</span>
            </a>
          </li>
        </ul>
      </li>
    </ul>
  </nav>
</template>

<script>
export default {
  name: 'OnThisPage',
  props: {
    containerSelector: {
      type: String,
      default: '.markdown-body'
    },
    levels: {
      type: Array,
      default: () => [2, 3, 4, 5, 6]
    },
    offset: {
      type: Number,
      default: 8
    }
  },
  data() {
    return {
      toc: [],
      activeId: '',
      _otpObserver: null,
      _otpObserverTimer: null,
      _otpPoller: null
    }
  },
  mounted() {
    this.buildToc();
    this.bindScrollSpy();
    
    this.$nextTick(() => {
      this.setupContainerObserver();
    });
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.onScrollSpy);
    window.removeEventListener('resize', this.onScrollSpy);
    this.cleanupObservers();
  },
  methods: {
    // 新增：重置目录状态
    resetToc() {
      this.toc = [];
      this.activeId = '';
      this.cleanupObservers();
      this.$nextTick(() => this.setupContainerObserver());
    },
    
    // 刷新目录内容
    refreshToc() {
      this.buildToc();
      this.onScrollSpy();
    },
    
    // 清理所有观察者
    cleanupObservers() {
      if (this._otpObserver) {
        this._otpObserver.disconnect();
        this._otpObserver = null;
      }
      if (this._otpObserverTimer) {
        clearTimeout(this._otpObserverTimer);
        this._otpObserverTimer = null;
      }
      if (this._otpPoller) {
        clearInterval(this._otpPoller);
        this._otpPoller = null;
      }
    },
    
    // 设置容器观察者
    setupContainerObserver() {
      // 多次重试确保目录能正确生成
      [0, 200, 500, 1000].forEach(delay => {
        setTimeout(() => this.refreshToc(), delay);
      });

      // 轮询检查容器是否存在
      const checkContainer = () => {
        const root = document.querySelector(this.containerSelector);
        if (root) {
          // 清除轮询
          if (this._otpPoller) {
            clearInterval(this._otpPoller);
            this._otpPoller = null;
          }

          // 监听容器变化
          if (!this._otpObserver) {
            this._otpObserver = new MutationObserver(() => {
              clearTimeout(this._otpObserverTimer);
              this._otpObserverTimer = setTimeout(() => this.refreshToc(), 100);
            });
            this._otpObserver.observe(root, {
              childList: true,
              subtree: true,
              attributes: true,
              attributeFilter: ['id']
            });
          }

          this.refreshToc();
        }
      };

      checkContainer();
      // 启动轮询（每200ms一次，直到找到容器）
      if (!this._otpPoller && !document.querySelector(this.containerSelector)) {
        this._otpPoller = setInterval(checkContainer, 200);
      }
    },
    
    getHeadingText(h) {
      try {
        const clone = h.cloneNode(true);
        clone.querySelectorAll('.heading-anchor')?.forEach(a => a.remove());
        let text = clone.textContent || '';
        text = text.replace(/\s*#\s*$/, '').trim();
        return text;
      } catch (e) {
        return (h.textContent || '').replace(/\s*#\s*$/, '').trim();
      }
    },
    
    buildToc() {
      const root = document.querySelector(this.containerSelector);
      if (!root) {
        this.toc = [];
        return;
      }
      
      const selector = this.levels.map(l => `h${l}`).join(',');
      const headings = Array.from(root.querySelectorAll(selector));
      
      if (headings.length === 0) {
        this.toc = [];
        return;
      }

      // 确保每个标题都有唯一id
      headings.forEach(h => {
        if (!h.id) {
          let safeId = h.textContent.trim().toLowerCase()
            .replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g, '')
            .replace(/\s+/g, '-');

          if (!safeId) {
            safeId = `section-${Math.random().toString(36).substring(2, 9)}`;
          }
          
          // 处理ID重复
          let finalId = safeId;
          let count = 1;
          while (document.getElementById(finalId)) {
            finalId = `${safeId}-${count++}`;
          }
          h.id = finalId;
        }
      });

      // 构建二级/三级嵌套目录
      const levelSet = new Set(this.levels);
      const topLevel = Math.min(...this.levels);
      const secondLevel = topLevel + 1;

      const toc = [];
      let currentTop = null;

      for (const h of headings) {
        const level = parseInt(h.tagName.substring(1), 10);
        if (!levelSet.has(level)) continue;

        const node = {
          id: h.id,
          text: this.getHeadingText(h),
          level,
          children: []
        };

        if (level === topLevel) {
          toc.push(node);
          currentTop = node;
        } else if (currentTop && level >= secondLevel) {
          currentTop.children.push(node);
        } else {
          toc.push(node);
        }
      }

      this.toc = toc;
    },
    
    bindScrollSpy() {
      window.addEventListener('scroll', this.onScrollSpy, { passive: true });
      window.addEventListener('resize', this.onScrollSpy);
      this.onScrollSpy();
    },
    
    onScrollSpy() {
      const root = document.querySelector(this.containerSelector);
      if (!root) return;
      
      const selector = this.levels.map(l => `h${l}`).join(',');
      const headings = Array.from(root.querySelectorAll(selector));
      if (headings.length === 0) {
        this.activeId = '';
        return;
      }
      
      const scrollY = window.scrollY || window.pageYOffset;

      let current = '';
      for (const h of headings) {
        const top = h.getBoundingClientRect().top + scrollY;
        if (top - this.offset <= scrollY + 1) {
          current = h.id;
        } else {
          break;
        }
      }
      this.activeId = current || (headings[0]?.id || '');
    },
    
    scrollToId(id) {
      const el = document.getElementById(id);
      if (!el) {
        return;
      }
      const top = el.getBoundingClientRect().top + window.scrollY - this.offset;
      window.scrollTo({
        top,
        behavior: 'smooth'
      });
    }
  }
}
</script>

<style scoped>
/* 样式保持不变 */
.on-this-page {
  --otp-border: var(--app-border);
  --otp-muted: var(--app-text-muted);
  --otp-bg-hover: var(--app-bg-tertiary);
  --otp-active-bg: var(--app-primary-bg-subtle);
  --otp-active: var(--app-primary);
  padding: 0.25rem 0;
  font-size: 0.95rem;
}

.otp-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  margin-bottom: 0.5rem;
}

.otp-title {
  font-weight: 600;
  color: var(--app-text-emphasis);
}

.otp-list {
  list-style: none;
  padding: 0;
  margin: 0;
}

.otp-item {
  margin: 0.125rem 0;
  padding-left: 0;
}

.otp-link {
  display: block;
  padding: 0.3rem 0.5rem;
  color: var(--app-text-muted);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: background-color 0.15s ease, color 0.15s ease;
}

.otp-link:hover {
  color: var(--app-primary);
}

.otp-item.active > .otp-link {
  color: var(--otp-active);
  font-weight: 600;
}

.otp-sublist {
  list-style: none;
  padding-left: 1.25rem;
  margin: 0.25rem 0 0.25rem 0;
}

.otp-subitem .otp-sublink {
  display: block;
  padding: 0.25rem 0.5rem;
  color: var(--app-text-muted);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: background-color 0.15s ease, color 0.15s ease;
  font-size: 0.9rem;
}

.otp-sublink:hover {
  color: var(--app-primary);
}

.otp-subitem.active > .otp-sublink {
  color: var(--otp-active);
  font-weight: 600;
}
</style>
