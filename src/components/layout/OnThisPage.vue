<!-- OnThisPage.vue (Unified: MarkdownToc + OnThisPage) -->
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
    // 正文容器选择器，默认存量 RenderMarkdown 输出使用 .markdown-body
    containerSelector: {
      type: String,
      default: '.markdown-body'
    },
    // 生成目录的标题层级
    levels: {
      type: Array,
      default: () => [2, 3, 4, 5, 6]
    },
    // 滚动偏移（如有固定头部）
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
      _otpObserverTimer: null
    }
  },
  mounted() {
    this.buildToc();
    this.bindScrollSpy();
    // 文章渲染后可能才有标题，延迟一次重建
    this.$nextTick(() => {
      setTimeout(() => {
        this.buildToc();
      }, 0);
      setTimeout(() => {
        this.buildToc();
      }, 200);
      setTimeout(() => {
        this.buildToc();
      }, 500);

      const root = document.querySelector(this.containerSelector);
      if (root && !this._otpObserver) {
        this._otpObserver = new MutationObserver(() => {
          clearTimeout(this._otpObserverTimer);
          this._otpObserverTimer = setTimeout(() => this.buildToc(), 100);
        });
        this._otpObserver.observe(root, { childList: true, subtree: true });
      }
    });
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.onScrollSpy);
    window.removeEventListener('resize', this.onScrollSpy);
    if (this._otpObserver) {
      this._otpObserver.disconnect();
      this._otpObserver = null;
    }
    if (this._otpObserverTimer) {
      clearTimeout(this._otpObserverTimer);
      this._otpObserverTimer = null;
    }
  },
  methods: {
    buildToc() {
      const root = document.querySelector(this.containerSelector);
      if (!root) {
        this.toc = [];
        return;
      }
      const selector = this.levels.map(l => `h${l}`).join(',');
      const headings = Array.from(root.querySelectorAll(selector));

      // 确保每个标题都有 id（若没有则根据文本生成）
      headings.forEach(h => {
        if (!h.id) {
          let safe = h.textContent.trim().toLowerCase()
            .replace(/[^\p{L}\p{N}\s-]/gu, '') // Keep Unicode letters, numbers, whitespace, hyphen
            .replace(/\s+/g, '-'); // Replace whitespace with hyphen

          if (!safe) {
            // Fallback for empty IDs (e.g., from headers with only symbols)
            safe = `section-${Math.random().toString(36).substring(2, 9)}`;
          }
          
          // 避免重复
          let id = safe;
          let n = 1;
          while (document.getElementById(id)) {
            id = `${safe}-${n++}`;
          }
          h.id = id;
        }
      });

      // 构建二级/三级嵌套（h2 -> children: h3/h4...）
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
          text: h.textContent.trim(),
          level,
          children: []
        };

        if (level === topLevel) {
          toc.push(node);
          currentTop = node;
        } else if (currentTop && level >= secondLevel) {
          currentTop.children.push(node);
        } else {
          // 如果没有 topLevel，直接平铺为顶层
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
        return;
      }
      const scrollY = window.scrollY || window.pageYOffset;

      let current = '';
      for (const h of headings) {
        const top = h.getBoundingClientRect().top + scrollY;
        // Add a 1px tolerance to handle rounding errors in scroll position
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
.on-this-page {
  --otp-border: var(--bs-border-color, #dee2e6);
  --otp-muted: var(--bs-secondary-color, #6c757d);
  --otp-bg-hover: var(--bs-tertiary-bg, #f8f9fa);
  --otp-active-bg: var(--bs-primary-bg-subtle, #e7f1ff);
  --otp-active: var(--bs-primary, #0d6efd);
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
  color: var(--bs-body-color);
}

.otp-list {
  list-style: none;
  padding: 0;
  margin: 0;
}

.otp-item {
  margin: 0.125rem 0;
  padding-left: 0.5rem;
}

.otp-link {
  display: block;
  padding: 0.3rem 0.5rem;
  color: var(--bs-body-color);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: background-color 0.15s ease, color 0.15s ease;
}

.otp-link:hover {
  color: var(--bs-primary);
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
  color: var(--bs-body-color);
  text-decoration: none;
  border-radius: 0.25rem;
  transition: background-color 0.15s ease, color 0.15s ease;
  font-size: 0.9rem;
}

.otp-sublink:hover {
  color: var(--bs-primary);
}

.otp-subitem.active > .otp-sublink {
  color: var(--otp-active);
  font-weight: 600;
}
</style>